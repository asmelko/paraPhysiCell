#include "random.h"

#include <random>

using namespace biofvm;

std::mt19937 generator;

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

void physicell::random::set_seed(unsigned int seed) { generator.seed(seed); }
