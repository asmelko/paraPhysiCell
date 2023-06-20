#include "random.h"

#include <random>

using namespace biofvm;

physicell::random& physicell::random::instance()
{
	static random instance;
	return instance;
}

real_t physicell::random::uniform(const real_t min, const real_t max)
{
	static std::mt19937 generator;
	std::uniform_real_distribution<real_t> distribution(min, max);
	return distribution(generator);
}
