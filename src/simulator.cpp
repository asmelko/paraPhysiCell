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

#define make_duration_pair(D) std::make_pair(D, #D)

using namespace biofvm;
using namespace physicell;

void simulator_durations::print_durations()
{
	std::array<std::pair<std::size_t, const char*>, 20> durations { make_duration_pair(diffusion),
																	make_duration_pair(secretion),
																	make_duration_pair(gradient),
																	make_duration_pair(custom_rules),
																	make_duration_pair(custom_interactions),
																	make_duration_pair(forces),
																	make_duration_pair(motility),
																	make_duration_pair(membrane),
																	make_duration_pair(spring),
																	make_duration_pair(position),
																	make_duration_pair(cell_cell),
																	make_duration_pair(container_mech),
																	make_duration_pair(mesh_mech),
																	make_duration_pair(neighbors_mech),
																	make_duration_pair(advance_phe),
																	make_duration_pair(container_phe),
																	make_duration_pair(mesh_phe),
																	make_duration_pair(neighbors_phe),
																	make_duration_pair(svg_save),
																	make_duration_pair(full_save) };

	std::size_t diff_total = 0;
	std::size_t mech_total = 0;
	std::size_t phen_total = 0;
	std::size_t save_total = 0;
	std::size_t total = 0;

	{
		index_t i = 0;
		for (; i < 2; i++)
			diff_total += durations[i].first;

		for (; i < 14; i++)
			mech_total += durations[i].first;

		for (; i < 18; i++)
			phen_total += durations[i].first;

		for (; i < 20; i++)
			save_total += durations[i].first;

		total = diff_total + mech_total + phen_total + save_total;
	}

	std::cout << "Duration ratios (total: " << (double)total / 1000 << "ms): "
			  << "Diffusion: " << (int)(diff_total / (double)total * 100)
			  << "% Mechanics: " << (int)(mech_total / (double)total * 100)
			  << "% Phenotype: " << (int)(phen_total / (double)total * 100)
			  << "% Save: " << (int)(save_total / (double)total * 100) << "%" << std::endl;

	std::sort(durations.begin(), durations.end(), [](auto& a, auto& b) { return a.first > b.first; });

	std::cout << "  ";
	for (index_t i = 0; i < 5; i++)
	{
		std::cout << durations[i].second << ": " << durations[i].first / 1000. << "ms ";
	}
	std::cout << std::endl;

	// zero them out
	diffusion = 0;
	secretion = 0;
	gradient = 0;
	custom_rules = 0;
	custom_interactions = 0;
	forces = 0;
	motility = 0;
	membrane = 0;
	spring = 0;
	position = 0;
	cell_cell = 0;
	container_mech = 0;
	mesh_mech = 0;
	neighbors_mech = 0;
	advance_phe = 0;
	container_phe = 0;
	mesh_phe = 0;
	neighbors_phe = 0;
	svg_save = 0;
	full_save = 0;
}

void simulator::initialize(environment& e, PhysiCell_Settings& settings)
{
	diffusion_solver_.initialize(e.m);
	mechanics_solver_.initialize(e);

	mechanics_step_interval_ = (index_t)std::round(e.mechanics_time_step / e.m.diffusion_time_step);
	phenotype_step_interval_ = (index_t)std::round(e.phenotype_time_step / e.m.diffusion_time_step);
	full_save_interval_ = (index_t)std::round(settings.full_save_interval / e.m.diffusion_time_step);
	svg_save_interval_ = (index_t)std::round(settings.SVG_save_interval / e.m.diffusion_time_step);

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

void simulator::simulate_diffusion_and_mechanics(environment& e, simulator_durations& durations,
												 index_t simulation_step, bool recompute_secretion_and_uptake)
{
	{
		// Compute diffusion:
		measure(diffusion_solver_.diffusion.solve(e.m), durations.diffusion);

		// Compute secretion and uptake:
		measure(diffusion_solver_.cell.simulate_secretion_and_uptake(e.m, recompute_secretion_and_uptake),
				durations.secretion);

		recompute_secretion_and_uptake = false;
	}

	if (simulation_step % mechanics_step_interval_ == 0)
	{
		// Compute gradient:
		measure(diffusion_solver_.gradient.solve(e.m), durations.gradient);

		// custom cell rules:
		measure(custom_cell_rules(e), durations.custom_rules);

		// Compute velocities and update the positions:
		{
			// custom attached cells adhesion:
			measure(evaluate_interactions(e), durations.custom_interactions);

			measure(mechanics_solver_.position.update_cell_forces(e), durations.forces);
			measure(mechanics_solver_.position.update_motility(e), durations.motility);
			measure(mechanics_solver_.position.update_basement_membrane_interactions(e), durations.membrane);
			measure(mechanics_solver_.position.update_spring_attachments(e), durations.spring);
			measure(mechanics_solver_.position.update_positions(e), durations.position);
		}

		// Standard cell-cell interactions:
		measure(mechanics_solver_.interactions.update_cell_cell_interactions(e), durations.cell_cell);

		// housekeeping
		{
			// Removal of flagged cells from data structures:
			measure(mechanics_solver_.containers.update_cell_container_for_mechanics(e), durations.container_mech);

			// Update mechanics mesh with new cell positions:
			measure(mechanics_solver_.containers.update_mechanics_mesh(e), durations.mesh_mech);

			// Update cells neighbors:
			measure(mechanics_solver_.position.update_cell_neighbors(e), durations.neighbors_mech);

			recompute_secretion_and_uptake = true;
		}
	}

	if (simulation_step % phenotype_step_interval_ == 0)
	{
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
}

void simulator::save_full(environment& e, const PhysiCell_Settings& settings, index_t simulation_step)
{
	if (settings.enable_full_saves == true)
	{
		std::stringstream ss;
		ss << settings.folder << "/output" << std::setfill('0') << std::setw(8)
		   << simulation_step / full_save_interval_;

		save_PhysiCell_to_MultiCellDS_v2(ss.str(), e);
	}
}

void simulator::save_svg(environment& e, const PhysiCell_Settings& settings,
						 const cell_coloring_funct_t& cell_coloring_function,
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

void simulator::run(environment& e, PhysiCell_Settings& settings, cell_coloring_funct_t cell_coloring_function,
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

		while (e.current_time < settings.max_time + 0.1 * e.m.diffusion_time_step)
		{
#pragma omp master
			{
				// save SVG plot if it's time
				if (simulation_step % svg_save_interval_ == 0)
				{
					measure(save_svg(e, settings, cell_coloring_function, substrate_coloring_function, simulation_step),
							durations.svg_save);
				}

				// save data if it's time.
				if (simulation_step % full_save_interval_ == 0)
				{
					measure(save_full(e, settings, simulation_step), durations.full_save);

					durations.print_durations();
					display_simulation_status(std::cout, e, settings);
				}
			}

			simulate_diffusion_and_mechanics(e, durations, simulation_step, recompute_secretion_and_uptake);

#pragma omp master
			e.current_time += e.m.diffusion_time_step;

			++simulation_step;
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
