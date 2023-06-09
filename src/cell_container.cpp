#include "cell_container.h"

#include <algorithm>

using namespace biofvm;
using namespace physicell;

cell_container::cell_container(microenvironment& m) : agent_container_templated<cell, cell_data>(m) {}

biofvm::agent_data& cell_container::get_agent_data() { return data_.agent_data; }

cell* cell_container::add_cell()
{
	agent_container_templated<cell, cell_data>::add_agent();

	return agents_.back().get();
}
