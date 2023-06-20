#pragma once

#include <BioFVM/agent_container_base.h>

#include "cell.h"
#include "cell_data.h"

namespace physicell {

class common_solver;
struct environment;

class cell_container_base
{
	friend common_solver;

protected:
	virtual cell_data& get_cell_data() = 0;
};

class cell_container : public biofvm::agent_container_common<cell, cell_data>, public cell_container_base
{
	virtual biofvm::agent_data& get_agent_data() override;
	virtual cell_data& get_cell_data() override;

public:
	cell_container(environment& e);
};

} // namespace physicell
