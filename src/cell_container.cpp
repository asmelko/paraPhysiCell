#include "cell_container.h"

#include <algorithm>

using namespace biofvm;
using namespace physicell;

cell_container::cell_container(environment& e) : agent_container_common<cell, cell_data>(e) {}

biofvm::agent_data& cell_container::get_agent_data() { return data_.agent_data; }

cell_data& cell_container::get_cell_data() { return data_; }

cell* cell_container::create_cell(cell_definition& definition)
{
	auto c = create();

	c->set_default(definition);

    return c;
}
