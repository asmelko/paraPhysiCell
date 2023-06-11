#include <chrono>
#include <iostream>

#include "BioFVM/solver.h"
#include "environment.h"
#include "solver/host/velocity_solver.h"

using namespace biofvm;
using namespace physicell;

void make_agents(microenvironment& m, index_t count, bool conflict)
{
	index_t x = 0, y = 0, z = 0;

	for (index_t i = 0; i < count; i++)
	{
		auto a = m.agents->create_agent();
		a->position()[0] = x;
		a->position()[1] = y;
		a->position()[2] = z;

		x += 20;
		if (x >= m.mesh.bounding_box_maxs[0])
		{
			x -= m.mesh.bounding_box_maxs[0];
			y += 20;
		}
		if (y >= m.mesh.bounding_box_maxs[1])
		{
			y -= m.mesh.bounding_box_maxs[1];
			z += 20;
		}
	}

	if (conflict)
	{
		auto a = m.agents->create_agent();
		a->position()[0] = 0;
		a->position()[1] = 0;
		a->position()[2] = 0;
	}
}

int main()
{
	cartesian_mesh mesh(3, { 0, 0, 0 }, { 5000, 5000, 5000 }, { 20, 20, 20 });

	real_t diffusion_time_step = 5;
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

	environment e(m);
	e.cell_definitions_count = cell_defs_count;
	m.agents = std::make_unique<cell_container>(e);

	make_agents(m, 2'000'000, true);

	solver s;

	s.initialize(m);

	for (index_t i = 0; i < 100; ++i)
	{
		std::size_t diffusion_duration, gradient_duration, secretion_duration, velocity_duration;
		{
			auto start = std::chrono::high_resolution_clock::now();

			s.diffusion.solve(m);

			auto end = std::chrono::high_resolution_clock::now();

			diffusion_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}

		{
			auto start = std::chrono::high_resolution_clock::now();

			s.gradient.solve(m);

			auto end = std::chrono::high_resolution_clock::now();

			gradient_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}

		{
			auto start = std::chrono::high_resolution_clock::now();

			s.cell.simulate_secretion_and_uptake(m, i % 10 == 0);

			auto end = std::chrono::high_resolution_clock::now();

			secretion_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}

		{
			auto start = std::chrono::high_resolution_clock::now();

			velocity_solver::solve(e);

			auto end = std::chrono::high_resolution_clock::now();

			velocity_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}

		std::cout << "Diffusion time: " << diffusion_duration << " ms,\t Gradient time: " << gradient_duration
				  << " ms,\t Secretion time: " << secretion_duration << " ms,\t Velocity time: " << velocity_duration
				  << " ms" << std::endl;
	}

	for (int i = 0; i < 5; i++)
	{
		s.diffusion.solve(m);
		s.gradient.solve(m);
	}
}
