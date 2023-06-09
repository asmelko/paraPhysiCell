#pragma once

#include <vector>

#include <BioFVM/agent_data.h>

namespace biofvm {
struct microenvironment;
}

namespace physicell {

struct cell_data
{
public:
	biofvm::agent_data agent_data;

	std::vector<biofvm::real_t> velocities;

	biofvm::index_t& agents_count; // references agent_data.agents_count

	biofvm::microenvironment& m;

	cell_data(biofvm::microenvironment& m);

	void add();
	void remove(biofvm::index_t index);
};

} // namespace physicell
