#pragma once

#include <BioFVM/types.h>

namespace physicell {

class random
{
public:
	static random& instance();

	biofvm::real_t uniform(const biofvm::real_t min = 0, const biofvm::real_t max = 1);

	void set_seed(unsigned int seed);
};

} // namespace physicell
