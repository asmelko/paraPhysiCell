#include "simulator.h"

#include <cmath>

using namespace biofvm;
using namespace physicell;

void simulator::initialize(environment& e)
{
	diffusion_solver_.initialize(e.m);
	mechanics_solver_.initialize(e);

	simulation_step_ = 0;
	mechanics_step_interval_ = (index_t)std::round(e.mechanics_time_step / e.m.time_step);
	recompute_secretion_and_uptake_ = false;
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
		mechanics_solver_.mechanics.update_mechanics_mesh(e);

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
		mechanics_solver_.mechanics.update_cell_container(e);

		recompute_secretion_and_uptake_ = true;
	}
}
