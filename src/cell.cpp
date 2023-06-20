#include "cell.h"

#include <BioFVM/microenvironment.h>

#include "cell_container.h"
#include "cell_data.h"
#include "environment.h"

using namespace biofvm;
using namespace physicell;

cell::cell(agent_id_t id, cell_data& data, index_t index)
	: agent(id, data.agent_data, index), data_(data), phenotype(data_, index)
{}

real_t* cell::velocity() { return data_.velocities.data() + index_ * data_.e.m.mesh.dims; }

index_t& cell::cell_definition_index() { return data_.cell_definition_indices[index_]; }

real_t& cell::simple_pressure() { return data_.simple_pressures[index_]; }

uint8_t& cell::is_movable() { return data_.is_movable[index_]; }

std::vector<cell*>& cell::neighbors() { return neighbors_; }
