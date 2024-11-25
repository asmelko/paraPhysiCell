#include "simulator.h"

#include <cmath>
#include <iomanip>
#include <sstream>

#include "original/core/standard_models.h"
#include "original/modules/BioFVM_MultiCellDS.h"
#include "original/modules/PhysiCell_MultiCellDS.h"
#include "original/modules/pathology.h"
#include "original/modules/timer.h"
#include "original/modules/various_output.h"

#define measure(F, D)                                                                                                  \
	{                                                                                                                  \
		auto start = std::chrono::high_resolution_clock::now();                                                        \
		F;                                                                                                             \
		auto end = std::chrono::high_resolution_clock::now();                                                          \
		D += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();                               \
	}

#define measure2(F, D)                                                                                                 \
	{                                                                                                                  \
		auto start = std::chrono::high_resolution_clock::now();                                                        \
		F;                                                                                                             \
		diffusion_solver_.wait_for_all();                                                                              \
		auto end = std::chrono::high_resolution_clock::now();                                                          \
		D += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();                               \
	}

using namespace biofvm;
using namespace physicell;

simulator::simulator(environment& e) : e(e) {}

void simulator::initialize(PhysiCell_Settings& settings)
{
	diffusion_solver_.initialize(e.m);
	mechanics_solver_.initialize(e);

	mechanics_step_interval_ = (index_t)std::round(e.mechanics_time_step / e.m.diffusion_time_step);
	phenotype_step_interval_ = (index_t)std::round(e.phenotype_time_step / e.m.diffusion_time_step);
	full_save_interval_ = (index_t)std::round(settings.full_save_interval / e.m.diffusion_time_step);
	svg_save_interval_ = (index_t)std::round(settings.SVG_save_interval / e.m.diffusion_time_step);
	max_time_ = (index_t)std::round(settings.max_time / e.m.diffusion_time_step);

	mechanics_solver_.containers.update_mechanics_mesh(e);
	mechanics_solver_.position.update_cell_neighbors(e);
}

void custom_cell_rules(environment& e)
{
	auto& cells = e.get_container();

#pragma omp for
	for (auto& cell : cells.agents())
	{
		if (cell->functions.custom_cell_rule)
		{
			cell->functions.custom_cell_rule(*cell);
		}
	}
}

void evaluate_interactions(environment& e)
{
	auto& cells = e.get_container();

#pragma omp for
	for (auto& cell : cells.agents())
	{
		if (cell->functions.contact_function)
		{
			for (auto& other : cell->state.attached_cells())
			{
				auto& other_cell = *cells.agents()[other];

				cell->functions.contact_function(*cell, other_cell);
			}
		}
	}
}

void advance_bundled_phenotype_functions(environment& e)
{
	auto& cells = e.get_container();

#pragma omp for
	for (auto& cell : cells.agents())
	{
		advance_bundled_phenotype_functions(*cell, e);
	}
}

void simulator::sync_data_host(simulator_durations& durations, biofvm::solvers::data_residency& residency)
{
	if constexpr (biofvm::solver::is_device_solver)
	{
		if (residency != biofvm::solvers::data_residency::host)
		{
#pragma omp master
			{
				measure(diffusion_solver_.load_data_from_solver(e.m), durations.host_sync);
			}
			residency = biofvm::solvers::data_residency::host;
#pragma omp barrier
		}
	}
}

void simulator::sync_data_device(simulator_durations& durations, biofvm::solvers::data_residency& residency)
{
	if (residency != biofvm::solvers::data_residency::device)
	{
		measure(diffusion_solver_.store_data_to_solver(e.m), durations.device_sync);
		residency = biofvm::solvers::data_residency::device;
	}
}

void simulator::simulate_diffusion_device(simulator_durations& durations, bool& recompute_secretion_and_uptake,
										  biofvm::solvers::data_residency& residency)
{
#pragma omp master
	{
		sync_data_device(durations, residency);

		// Compute diffusion:
		measure(diffusion_solver_.diffusion.solve(e.m), durations.diffusion);

		// Compute secretion and uptake:
		measure(diffusion_solver_.cell.simulate_secretion_and_uptake(e.m, recompute_secretion_and_uptake),
				durations.secretion);
	}
	residency = biofvm::solvers::data_residency::device;
}

void simulator::simulate_diffusion(simulator_durations& durations, bool& recompute_secretion_and_uptake,
								   biofvm::solvers::data_residency& residency)
{
	if constexpr (biofvm::solver::is_device_solver)
	{
		simulate_diffusion_device(durations, recompute_secretion_and_uptake, residency);
	}
	else
	{
		// Compute diffusion:
		measure(diffusion_solver_.diffusion.solve(e.m), durations.diffusion);

		// Compute secretion and uptake:
		measure(diffusion_solver_.cell.simulate_secretion_and_uptake(e.m, recompute_secretion_and_uptake),
				durations.secretion);
	}

	recompute_secretion_and_uptake = false;
}

void simulator::simulate_mechanics(simulator_durations& durations, bool& recompute_secretion_and_uptake,
								   biofvm::solvers::data_residency& residency)
{
	sync_data_host(durations, residency);

	// Compute gradient:
	// measure(diffusion_solver_.gradient.solve(e.m), durations.gradient);

	// custom cell rules:
	// measure(custom_cell_rules(e), durations.custom_rules);

	// Compute velocities and update the positions:
	{
		// custom attached cells adhesion:
		// measure(evaluate_interactions(e), durations.custom_interactions);

		measure(mechanics_solver_.position.update_cell_forces_new(e), durations.forces);
		measure(mechanics_solver_.position.update_motility(e), durations.motility);
		measure(mechanics_solver_.position.update_basement_membrane_interactions(e), durations.membrane);
		// measure(mechanics_solver_.position.update_spring_attachments(e), durations.spring);
		measure(mechanics_solver_.position.update_positions_new(e), durations.position);
	}

	// Standard cell-cell interactions:
	// measure(mechanics_solver_.interactions.update_cell_cell_interactions(e), durations.cell_cell);

	// housekeeping
	{
		// Removal of flagged cells from data structures:
		// measure(mechanics_solver_.containers.update_cell_container_for_mechanics(e), durations.container_mech);

		// Update mechanics mesh with new cell positions:
		// measure(mechanics_solver_.containers.update_mechanics_mesh(e), durations.mesh_mech);

		// Update cells neighbors:
		// measure(mechanics_solver_.position.update_cell_neighbors(e), durations.neighbors_mech);

		recompute_secretion_and_uptake = true;
	}
}

void simulator::simulate_phenotype(simulator_durations& durations, bool& recompute_secretion_and_uptake,
								   biofvm::solvers::data_residency& residency)
{
	sync_data_host(durations, residency);

	// Update phenotype:
	measure(::advance_bundled_phenotype_functions(e), durations.advance_phe);

	// housekeeping
	{
		// Removal and division of flagged cells in data structures:
		measure(mechanics_solver_.containers.update_cell_container_for_phenotype(e, diffusion_solver_.cell),
				durations.container_phe);

		// Update mechanics mesh with new cell positions:
		measure(mechanics_solver_.containers.update_mechanics_mesh(e), durations.mesh_phe);

		// Update cells neighbors:
		measure(mechanics_solver_.position.update_cell_neighbors(e), durations.neighbors_phe);

		recompute_secretion_and_uptake = true;
	}
}

void simulator::save(simulator_durations& durations, PhysiCell_Settings& settings,
					 const cell_coloring_funct_t& cell_coloring_function,
					 const substrate_coloring_funct_t& substrate_coloring_function, index_t simulation_step)
{
#pragma omp master
	{
		// save SVG plot if it's time
		if (simulation_step % svg_save_interval_ == 0)
		{
			measure(save_svg(settings, cell_coloring_function, substrate_coloring_function, simulation_step),
					durations.svg_save);
		}

		// save data if it's time.
		if (simulation_step % full_save_interval_ == 0)
		{
			measure(save_full(settings, simulation_step), durations.full_save);

			durations.print_durations();
			display_simulation_status(std::cout, e, settings);
		}
	}
}

void simulator::save_full(const PhysiCell_Settings& settings, index_t simulation_step)
{
	if (settings.enable_full_saves == true)
	{
		std::stringstream ss;
		ss << settings.folder << "/output" << std::setfill('0') << std::setw(8)
		   << simulation_step / full_save_interval_;

		save_PhysiCell_to_MultiCellDS_v2(ss.str(), e);
	}
}

void simulator::save_svg(const PhysiCell_Settings& settings, const cell_coloring_funct_t& cell_coloring_function,
						 const substrate_coloring_funct_t& substrate_coloring_function, index_t simulation_step)
{
	if (settings.enable_SVG_saves == true)
	{
		std::stringstream ss;
		ss << settings.folder << "/snapshot" << std::setfill('0') << std::setw(8)
		   << simulation_step / svg_save_interval_ << ".svg";

		SVG_plot(ss.str(), e, settings, 0.0, e.current_time, cell_coloring_function, substrate_coloring_function);
	}
}

void simulator::run(PhysiCell_Settings& settings, cell_coloring_funct_t cell_coloring_function,
					substrate_coloring_funct_t substrate_coloring_function)
{
	// set MultiCellDS save options

	set_save_biofvm_mesh_as_matlab(true);
	set_save_biofvm_data_as_matlab(true);
	set_save_biofvm_cell_data(true);
	set_save_biofvm_cell_data_as_custom_matlab(true);

	// save a simulation snapshot

	save_PhysiCell_to_MultiCellDS_v2(settings.folder + "/initial", e);

	// save a quick SVG cross section through z = 0, after setting its
	// length bar to 200 microns

	PhysiCell_SVG_options.length_bar = 200;

	// for simplicity, set a pathology coloring function

	SVG_plot(settings.folder + "/initial.svg", e, settings, 0.0, e.current_time, cell_coloring_function,
			 substrate_coloring_function);

	create_plot_legend(settings.folder + "/legend.svg", cell_coloring_function, e);

	// set the performance timers

	RUNTIME_TIC();
	TIC();

#pragma omp parallel
	{
		index_t simulation_step = 0;
		simulator_durations durations;
		bool recompute_secretion_and_uptake = true;
		biofvm::solvers::data_residency data_residency = biofvm::solvers::data_residency::host;

		while (simulation_step <= max_time_)
		{
			save(durations, settings, cell_coloring_function, substrate_coloring_function, simulation_step);

			// called each time because one simulation step is equal to one diffusion time step
			{
				// simulate_diffusion(durations, recompute_secretion_and_uptake, data_residency);
			}

			if (simulation_step % mechanics_step_interval_ == 0)
			{
				simulate_mechanics(durations, recompute_secretion_and_uptake, data_residency);
			}

			if (simulation_step % phenotype_step_interval_ == 0)
			{
				// simulate_phenotype(durations, recompute_secretion_and_uptake, data_residency);
			}

			// run custom simulations
			for (const auto& [interval, custom_simulate] : custom_simulations_)
			{
				if (simulation_step % interval == 0)
				{
					custom_simulate(e, durations, recompute_secretion_and_uptake, data_residency);
				}
			}

#pragma omp master
			e.current_time += e.m.diffusion_time_step;

			++simulation_step;
			// std::exit(1);
		}
	}

	// save a final simulation snapshot

	save_PhysiCell_to_MultiCellDS_v2(settings.folder + "/final", e);

	SVG_plot(settings.folder + "/final.svg", e, settings, 0.0, e.current_time, cell_coloring_function,
			 substrate_coloring_function);

	// timer

	std::cout << std::endl << "Total simulation runtime: " << std::endl;
	display_stopwatch_value(std::cout, runtime_stopwatch_value());
}

void simulator::add_simulation_step(biofvm::real_t time_step, simulate_func_t&& simulate_f)
{
	index_t step_interval = (index_t)std::round(time_step / e.m.diffusion_time_step);

	custom_simulations_.emplace_back(step_interval, [=, this](environment& e, simulator_durations& durations,
															  bool& recompute_secretion_and_uptake,
															  biofvm::solvers::data_residency& residency) {
		recompute_secretion_and_uptake = true;
		sync_data_host(durations, residency);
		simulate_f(e);
	});
}
