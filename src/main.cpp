#include <chrono>
#include <iostream>

#include <BioFVM/solver.h>

using namespace biofvm;

int main()
{
	cartesian_mesh mesh(3, { 0, 0, 0 }, { 500, 500, 500 }, { 20, 20, 20 });

	real_t diffusion_time_step = 5;
	index_t substrates_count = 4;

	auto diff_coefs = std::make_unique<real_t[]>(substrates_count);
	diff_coefs[0] = 4;
	auto decay_rates = std::make_unique<real_t[]>(substrates_count);
	decay_rates[0] = 5;

	auto initial_conds = std::make_unique<real_t[]>(substrates_count);
	initial_conds[0] = 0;

	microenvironment m(mesh, substrates_count, diffusion_time_step, initial_conds.get());

	m.diffustion_coefficients = std::move(diff_coefs);
	m.decay_rates = std::move(decay_rates);

	solver s;

	s.initialize(m);

	for (index_t i = 0; i < 100; ++i)
	{
		std::size_t diffusion_duration, gradient_duration;
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

		std::cout << "Diffusion time: " << diffusion_duration << " ms,\t Gradient time: " << gradient_duration << " ms"
				  << std::endl;
	}

	for (int i = 0; i < 5; i++)
	{
		s.diffusion.solve(m);
		s.gradient.solve(m);
	}
}
