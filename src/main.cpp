#include <chrono>
#include <iostream>

#include "BioFVM/solver.h"
#include "environment.h"
#include "solver/host/mechanics_solver.h"
#include "solver/host/position_solver.h"

using namespace biofvm;
using namespace physicell;

#define measure(F, D)                                                                                                  \
	{                                                                                                                  \
		auto start = std::chrono::high_resolution_clock::now();                                                        \
		F;                                                                                                             \
		auto end = std::chrono::high_resolution_clock::now();                                                          \
		D = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();                                \
	}

void make_agents(environment& e, index_t count, bool conflict)
{
	index_t x = 0, y = 0, z = 0;

	for (index_t i = 0; i < count; i++)
	{
		auto a = e.cast_container<cell_container>().create();
		a->position()[0] = x + 10;
		a->position()[1] = y + 10;
		a->position()[2] = z + 10;

		a->phenotype.geometry.radius() = 15;
		a->phenotype.mechanics.relative_maximum_adhesion_distance() = 1;

		a->is_movable() = 1;
		a->phenotype.motility.is_motile() = 1;

		a->phenotype.mechanics.attachment_rate() = 100;
		a->phenotype.mechanics.detachment_rate() = 100;
		a->phenotype.mechanics.maximum_number_of_attachments() = 12;
		a->phenotype.mechanics.cell_adhesion_affinities()[0] = 10;

		a->phenotype.mechanics.cell_BM_repulsion_strength() = 100;

		x += e.m.mesh.voxel_shape[0];
		if (x >= e.m.mesh.bounding_box_maxs[0])
		{
			x -= e.m.mesh.bounding_box_maxs[0];
			y += e.m.mesh.voxel_shape[1];
		}
		if (y >= e.m.mesh.bounding_box_maxs[1])
		{
			y -= e.m.mesh.bounding_box_maxs[1];
			z += e.m.mesh.voxel_shape[2];
		}
	}

	if (conflict)
	{
		auto a = e.m.agents->create_agent();
		a->position()[0] = 10;
		a->position()[1] = 10;
		a->position()[2] = 10;
	}
}

int main()
{
	cartesian_mesh mesh(3, { 0, 0, 0 }, { 5000, 5000, 5000 }, { 20, 20, 20 });
	cartesian_mesh mechanics_mesh(3, { 0, 0, 0 }, { 5000, 5000, 5000 }, { 40, 40, 40 });

	real_t diffusion_time_step = 0.01;
	index_t substrates_count = 4;
	index_t cell_defs_count = 3;

	auto diff_coefs = std::make_unique<real_t[]>(substrates_count);
	diff_coefs[0] = 4;
	auto decay_rates = std::make_unique<real_t[]>(substrates_count);
	decay_rates[0] = 5;

	auto initial_conds = std::make_unique<real_t[]>(substrates_count);
	initial_conds[0] = 0;

	microenvironment m(mesh, substrates_count, diffusion_time_step, initial_conds.get());

	m.diffustion_coefficients = std::move(diff_coefs);
	m.decay_rates = std::move(decay_rates);
	m.compute_internalized_substrates = true;

	environment e(m, mechanics_mesh);
	e.cell_definitions_count = cell_defs_count;
	m.agents = std::make_unique<cell_container>(e);
	e.mechanics_time_step = 0.1;

	make_agents(e, 2'000'000, true);

	solver s;

	s.initialize(m);

	if (true)
		for (index_t i = 0; i < 100; ++i)
		{
			std::size_t diffusion_duration, gradient_duration, secretion_duration, velocity_update_mesh,
				velocity_interactions_duration, velocity_motility_duration, velocity_membrane_duration,
				velocity_attachments_duration, position_duration;

			measure(s.diffusion.solve(m), diffusion_duration);
			measure(s.gradient.solve(m), gradient_duration);
			measure(s.cell.simulate_secretion_and_uptake(m, i % 10 == 0), secretion_duration);
			measure(mechanics_solver::update_mechanics_mesh(e), velocity_update_mesh);
			measure(position_solver::update_cell_velocities_and_neighbors(e), velocity_interactions_duration);
			measure(position_solver::update_motility(e), velocity_motility_duration);
			measure(position_solver::update_basement_membrane_interactions(e), velocity_membrane_duration);
			measure(position_solver::update_spring_attachments(e), velocity_attachments_duration);
			measure(position_solver::update_positions(e), position_duration);

			std::cout << "Diffusion time: " << diffusion_duration << " ms,\t Gradient time: " << gradient_duration
					  << " ms,\t Secretion time: " << secretion_duration
					  << " ms,\t Positions[mech,cells,motil,membr,spring,upd] time: " << velocity_update_mesh << ","
					  << velocity_interactions_duration << "," << velocity_motility_duration << ","
					  << velocity_membrane_duration << "," << velocity_attachments_duration << "," << position_duration
					  << " ms" << std::endl;
		}

	for (int i = 0; i < 5; i++)
	{
		mechanics_solver::update_mechanics_mesh(e);
		position_solver::update_cell_velocities_and_neighbors(e);
		position_solver::update_motility(e);
		position_solver::update_basement_membrane_interactions(e);
		position_solver::update_spring_attachments(e);
		position_solver::update_positions(e);
	}
}
