#pragma once

#include <BioFVM/types.h>

namespace physicell {

class random
{
public:
	static random& instance();

	biofvm::real_t uniform(const biofvm::real_t min = 0, const biofvm::real_t max = 1);

	biofvm::real_t normal(const biofvm::real_t mean = 0, const biofvm::real_t std = 1);

	void set_seed(unsigned int seed);

	void set_random_seed();
};

} // namespace physicell
