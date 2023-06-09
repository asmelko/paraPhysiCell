#include "cell_container.h"

#include <algorithm>

using namespace biofvm;
using namespace physicell;

cell_container::cell_container(microenvironment& m, index_t cell_definitions_count)
	: agent_container_common<cell, cell_data>(m, cell_definitions_count)
{}

biofvm::agent_data& cell_container::get_agent_data() { return data_.agent_data; }

cell* cell_container::add_cell()
{
	agent_container_common<cell, cell_data>::add_agent();

	return agents_.back().get();
}
