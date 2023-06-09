#include "cell_data.h"

#include <BioFVM/microenvironment.h>

using namespace biofvm;
using namespace physicell;

cell_data::cell_data(microenvironment& m) : agent_data(m), agents_count(agent_data.agents_count), m(m) {}

void cell_data::add()
{
	agent_data.add();

	for (index_t i = 0; i < m.mesh.dims; i++)
	{
		velocities.push_back(0);
	}
}

void cell_data::remove(index_t index)
{
	agent_data.remove(index);

	if (index == agents_count)
	{
		return;
	}

	for (index_t i = 0; i < m.mesh.dims; i++)
	{
		velocities[index * m.mesh.dims + i] = velocities[agents_count * m.mesh.dims + i];
	}
}
