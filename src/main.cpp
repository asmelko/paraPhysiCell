#include <chrono>
#include <iostream>

#include <BioFVM/solver.h>

#include "environment.h"
#include "original/core/standard_models.h"
#include "solver/host/containers_solver.h"
#include "solver/host/interactions_solver.h"
#include "solver/host/position_solver.h"

using namespace biofvm;
using namespace physicell;

#define measure(F, D)                                                                                                  \
	auto D##start = std::chrono::high_resolution_clock::now();                                                         \
	F;                                                                                                                 \
	auto D##end = std::chrono::high_resolution_clock::now();                                                           \
	D = std::chrono::duration_cast<std::chrono::milliseconds>(D##end - D##start).count();

void make_agents(environment& e, index_t count, bool conflict)
{
	index_t x = 0, y = 0, z = 0;

	for (index_t i = 0; i < count; i++)
	{
		auto a = e.get_container().create_cell(e.cell_defaults());
		a->position()[0] = x + 10;
		a->position()[1] = y + 10;
		a->position()[2] = z + 10;

		a->phenotype.geometry.radius() = 15;
		a->phenotype.mechanics.relative_maximum_adhesion_distance() = 1;
		a->phenotype.volume.nuclear_fluid() = 100;
		a->phenotype.volume.nuclear_solid() = 100;
		a->phenotype.volume.cytoplasmic_fluid() = 100;
		a->phenotype.volume.cytoplasmic_solid() = 100;
		a->phenotype.volume.fluid() = 200;
		a->phenotype.volume.solid() = 200;
		a->phenotype.volume.total() = 400;

		a->is_movable() = 1;
		a->phenotype.motility.is_motile() = 1;

		a->phenotype.mechanics.attachment_rate() = .1;
		a->phenotype.mechanics.detachment_rate() = .1;
		a->phenotype.mechanics.maximum_number_of_attachments() = 12;
		a->phenotype.mechanics.cell_adhesion_affinities()[0] = .1;

		a->phenotype.mechanics.cell_BM_repulsion_strength() = 100;

		a->phenotype.cell_interactions.attack_rates()[0] = .01;
		a->phenotype.cell_interactions.fusion_rates()[0] = .01;
		a->phenotype.cell_interactions.live_phagocytosis_rates()[0] = .01;

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
		auto a = e.get_container().create();
		a->position()[0] = 10;
		a->position()[1] = 10;
		a->position()[2] = 10;
		a->phenotype.volume.nuclear_fluid() = 100;
		a->phenotype.volume.nuclear_solid() = 100;
		a->phenotype.volume.cytoplasmic_fluid() = 100;
		a->phenotype.volume.cytoplasmic_solid() = 100;
		a->phenotype.volume.fluid() = 200;
		a->phenotype.volume.solid() = 200;
		a->phenotype.volume.total() = 400;
	}
}

int main()
{
	std::size_t microenv_init_duration, env_init_duration, cells_init_duration, solver_init_duration;

	cartesian_mesh mesh(3, { 0, 0, 0 }, { 10000, 10000, 10000 }, { 20, 20, 20 });

	real_t diffusion_time_step = 0.01;
	index_t substrates_count = 2;
	index_t cell_defs_count = 3;

	auto diff_coefs = std::make_unique<real_t[]>(substrates_count);
	diff_coefs[0] = 4;
	auto decay_rates = std::make_unique<real_t[]>(substrates_count);
	decay_rates[0] = 5;

	auto initial_conds = std::make_unique<real_t[]>(substrates_count);
	initial_conds[0] = 0;

	measure(microenvironment m(mesh, substrates_count, diffusion_time_step, initial_conds.get()),
			microenv_init_duration);
	m.diffusion_coefficients = std::move(diff_coefs);
	m.decay_rates = std::move(decay_rates);
	m.compute_internalized_substrates = true;

	measure(environment e(m, cell_defs_count, { 40, 40, 40 }), env_init_duration);
	e.cell_definitions_count = cell_defs_count;
	m.agents = std::make_unique<cell_container>(e);
	e.mechanics_time_step = 0.1;
	e.virtual_wall_at_domain_edges = true;
	e.automated_spring_adhesion = true;
	initialize_default_cell_definition(e);

	measure(make_agents(e, 2'000'000, true), cells_init_duration);

	solver s;
	interactions_solver is;
	containers_solver cs;
	cs.initialize(e);

	measure(s.initialize(m), solver_init_duration);

	std::cout << "Init time[m, e, c, s]: " << microenv_init_duration << "," << env_init_duration << ","
			  << cells_init_duration << "," << solver_init_duration << " ms" << std::endl;

	if (true)
		for (index_t i = 0; i < 100; ++i)
		{
			std::size_t diffusion_duration, gradient_duration, secretion_duration, velocity_update_mesh,
				velocity_forces_duration, neighbors_duration, velocity_motility_duration, velocity_membrane_duration,
				velocity_attachments_duration, position_duration, interactions_duration, delete_duration;

#pragma omp parallel
			{
				measure(s.diffusion.solve(m), diffusion_duration);
				measure(s.gradient.solve(m), gradient_duration);
				measure(s.cell.simulate_secretion_and_uptake(m, i % 10 == 0), secretion_duration);
				measure(cs.update_mechanics_mesh(e), velocity_update_mesh);
				measure(position_solver::update_cell_neighbors(e), neighbors_duration);
				measure(position_solver::update_cell_forces(e), velocity_forces_duration);
				measure(position_solver::update_motility(e), velocity_motility_duration);
				measure(position_solver::update_basement_membrane_interactions(e), velocity_membrane_duration);
				measure(position_solver::update_spring_attachments(e), velocity_attachments_duration);
				measure(position_solver::update_positions(e), position_duration);
				measure(is.update_cell_cell_interactions(e), interactions_duration);
				measure(cs.update_cell_container_for_mechanics(e), delete_duration);
			}

			std::cout << "Diffusion time: " << diffusion_duration << " ms,\t Gradient time: " << gradient_duration
					  << " ms,\t Secretion time: " << secretion_duration
					  << " ms,\t Positions[mech,cells,nei,motil,membr,spring,upd] time: " << velocity_update_mesh << ","
					  << velocity_forces_duration << "," << neighbors_duration << "," << velocity_motility_duration
					  << "," << velocity_membrane_duration << "," << velocity_attachments_duration << ","
					  << position_duration << " ms,\t Interactions time: " << interactions_duration
					  << " ms,\t Delete time: " << delete_duration << " ms" << std::endl;

			std::cout << "Number of cells: " << e.get_container().data().agents_count << std::endl;
		}

	for (index_t i = 0; i < 5; i++)
	{
		cs.update_mechanics_mesh(e);
		position_solver::update_cell_neighbors(e);
		position_solver::update_cell_forces(e);
		position_solver::update_motility(e);
		position_solver::update_basement_membrane_interactions(e);
		position_solver::update_spring_attachments(e);
		position_solver::update_positions(e);
	}
}
