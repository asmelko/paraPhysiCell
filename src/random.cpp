#include "random.h"

#include <omp.h>
#include <random>

using namespace biofvm;

thread_local std::mt19937 generator;

physicell::random& physicell::random::instance()
{
	static random instance;
	return instance;
}

real_t physicell::random::uniform(const real_t min, const real_t max)
{
	std::uniform_real_distribution<real_t> distribution(min, max);
	return distribution(generator);
}

real_t physicell::random::normal(const real_t mean, const real_t std)
{
	std::normal_distribution<real_t> distribution(mean, std);
	return distribution(generator);
}

void physicell::random::set_seed(unsigned int seed)
{
#ifdef _OPENMP
	std::vector<unsigned int> initial_sequence(omp_get_num_threads());

	for (int i = 0; i < omp_get_num_threads(); i++)
		initial_sequence[i] = seed + i;

	std::seed_seq seq(initial_sequence.begin(), initial_sequence.end());

	std::vector<unsigned int> seeds(omp_get_num_threads());
	seq.generate(seeds.begin(), seeds.end());

	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		generator.seed(seeds[id]);
	}
#else
	generator.seed(seed);
#endif
}

void physicell::random::set_random_seed()
{
	std::random_device rd;
	set_seed(rd());
}
