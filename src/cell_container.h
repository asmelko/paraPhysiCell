#pragma once

#include <BioFVM/agent_container_base.h>

#include "cell.h"
#include "cell_data.h"

namespace physicell {

class cell_container : public biofvm::agent_container_common<cell, cell_data>
{
	virtual biofvm::agent_data& get_agent_data() override;

public:
	cell_container(biofvm::microenvironment& m, biofvm::index_t cell_definitions_count);

	cell* add_cell();
};

} // namespace physicell
