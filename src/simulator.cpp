#include "simulator.h"

#include <cmath>

#include "environment.h"
#include "original/standard_models.h"

using namespace biofvm;
using namespace physicell;

void simulator::initialize(environment& e)
{
	diffusion_solver_.initialize(e.m);
	mechanics_solver_.initialize(e);

	simulation_step_ = 0;
	mechanics_step_interval_ = (index_t)std::round(e.mechanics_time_step / e.m.time_step);
	phenotype_step_interval_ = (index_t)std::round(e.phenotype_time_step / e.m.time_step);
	recompute_secretion_and_uptake_ = false;
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

		// Update mechanics mesh with new cell positions:
		mechanics_solver_.containers.update_mechanics_mesh(e);

		// custom attached cells adhesion:
		evaluate_interactions(e);

		// custom cell rules:
		custom_cell_rules(e);

		// Compute velocities and update the positions:
		{
			mechanics_solver_.position.update_cell_velocities_and_neighbors(e);
			mechanics_solver_.position.update_motility(e);
			mechanics_solver_.position.update_basement_membrane_interactions(e);
			mechanics_solver_.position.update_spring_attachments(e);
			mechanics_solver_.position.update_positions(e);
		}

		// Standard cell-cell interactions:
		mechanics_solver_.interactions.update_cell_cell_interactions(e);

		// Removal of flagged cells from data structures:
		mechanics_solver_.containers.update_cell_container_for_mechanics(e);

		recompute_secretion_and_uptake_ = true;
	}

	if (simulation_step_ % phenotype_step_interval_ == 0)
	{
		// Update phenotype:
		advance_bundled_phenotype_functions(e);

		// Removal and division of flagged cells in data structures:
		mechanics_solver_.containers.update_cell_container_for_phenotype(e, diffusion_solver_.cell);
	}
}
