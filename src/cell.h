#pragma once

#include <BioFVM/agent.h>

namespace physicell {

struct cell_data;
class cell_container;

class cell : public biofvm::agent
{
	friend cell_container;

	cell_data& data_;

public:
	cell(biofvm::agent_id_t id, cell_data& data, biofvm::index_t index);

	biofvm::real_t* velocities();
};

} // namespace physicell
