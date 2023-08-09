#include "simulator.h"

#include <cmath>
#include <iomanip>
#include <sstream>

#include "environment.h"
#include "original/modules/BioFVM_MultiCellDS.h"
#include "original/modules/PhysiCell_MultiCellDS.h"
#include "original/modules/pathology.h"
#include "original/modules/timer.h"
#include "original/modules/various_output.h"
#include "original/standard_models.h"

using namespace biofvm;
using namespace physicell;

void simulator::initialize(environment& e, PhysiCell_Settings& settings)
{
	diffusion_solver_.initialize(e.m);
	mechanics_solver_.initialize(e);

	simulation_step_ = 0;

	mechanics_step_interval_ = (index_t)std::round(e.mechanics_time_step / e.m.diffusion_time_step);
	phenotype_step_interval_ = (index_t)std::round(e.phenotype_time_step / e.m.diffusion_time_step);
	full_save_interval_ = (index_t)std::round(settings.full_save_interval / e.m.diffusion_time_step);
	svg_save_interval_ = (index_t)std::round(settings.SVG_save_interval / e.m.diffusion_time_step);

	recompute_secretion_and_uptake_ = true;
	
	mechanics_solver_.containers.update_cell_container_for_mechanics(e);
	mechanics_solver_.containers.update_mechanics_mesh(e);
}

void custom_cell_rules(environment& e)
{
	auto& cells = e.cast_container<cell_container>();

	for (auto& cell : cells.agents())
	{
		if (cell->functions.custom_cell_rule)
		{
			cell->functions.custom_cell_rule(*cell);
		}
	}
}

void simulator::simulate_diffusion_and_mechanics(environment& e)
{
	++simulation_step_;

	{
		// Compute diffusion:
		diffusion_solver_.diffusion.solve(e.m);

		// Compute secretion and uptake:
		diffusion_solver_.cell.simulate_secretion_and_uptake(e.m, recompute_secretion_and_uptake_);

		recompute_secretion_and_uptake_ = false;
	}

	if (simulation_step_ % mechanics_step_interval_ == 0)
	{
		// Compute gradient:
		diffusion_solver_.gradient.solve(e.m);

		// custom cell rules:
		custom_cell_rules(e);

		// Compute velocities and update the positions:
		{
			// custom attached cells adhesion:
			evaluate_interactions(e);

			mechanics_solver_.position.update_cell_forces(e);
			mechanics_solver_.position.update_motility(e);
			mechanics_solver_.position.update_basement_membrane_interactions(e);
			mechanics_solver_.position.update_spring_attachments(e);
			mechanics_solver_.position.update_positions(e);
		}

		// Standard cell-cell interactions:
		mechanics_solver_.interactions.update_cell_cell_interactions(e);

		// housekeeping
		{
			// Removal of flagged cells from data structures:
			mechanics_solver_.containers.update_cell_container_for_mechanics(e);

			// Update mechanics mesh with new cell positions:
			mechanics_solver_.containers.update_mechanics_mesh(e);

			// Update cells neighbors:
			mechanics_solver_.position.update_cell_neighbors(e);
		}

		recompute_secretion_and_uptake_ = true;
	}

	if (simulation_step_ % phenotype_step_interval_ == 0)
	{
		// Update phenotype:
		advance_bundled_phenotype_functions(e);

		// housekeeping
		{
			// Removal and division of flagged cells in data structures:
			mechanics_solver_.containers.update_cell_container_for_phenotype(e, diffusion_solver_.cell);

			// Update cells neighbors:
			mechanics_solver_.position.update_cell_neighbors(e);
		}
	}
}

void simulator::run(environment& e, PhysiCell_Settings& settings, cell_coloring_funct_t cell_coloring_function)
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

	SVG_plot(settings.folder + "/initial.svg", e, settings, 0.0, e.current_time, cell_coloring_function);

	create_plot_legend(settings.folder + "/legend.svg", cell_coloring_function, e);

	// set the performance timers

	RUNTIME_TIC();
	TIC();

	while (e.current_time < settings.max_time + 0.1 * e.m.diffusion_time_step)
	{
		// save data if it's time.
		if (simulation_step_ % full_save_interval_ == 0)
		{
			display_simulation_status(std::cout, e, settings);

			if (settings.enable_full_saves == true)
			{
				std::stringstream ss;
				ss << settings.folder << "/output" << std::setfill('0') << std::setw(8)
				   << simulation_step_ / full_save_interval_;

				save_PhysiCell_to_MultiCellDS_v2(ss.str(), e);
			}
		}

		// save SVG plot if it's time
		if (simulation_step_ % svg_save_interval_ == 0)
		{
			if (settings.enable_SVG_saves == true)
			{
				std::stringstream ss;
				ss << settings.folder << "/snapshot" << std::setfill('0') << std::setw(8)
				   << simulation_step_ / svg_save_interval_ << ".svg";

				SVG_plot(ss.str(), e, settings, 0.0, e.current_time, cell_coloring_function);
			}
		}

		simulate_diffusion_and_mechanics(e);

		e.current_time += e.m.diffusion_time_step;
	}

	// save a final simulation snapshot

	save_PhysiCell_to_MultiCellDS_v2(settings.folder + "/final", e);

	SVG_plot(settings.folder + "/final.svg", e, settings, 0.0, e.current_time, cell_coloring_function);

	// timer

	std::cout << std::endl << "Total simulation runtime: " << std::endl;
	display_stopwatch_value(std::cout, runtime_stopwatch_value());
}
