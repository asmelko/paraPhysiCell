#include "cell.h"

#include <BioFVM/microenvironment.h>

#include "cell_data.h"
#include "environment.h"

using namespace biofvm;
using namespace physicell;

cell_state_t::cell_state_t(cell_data& data, index_t index) : phenotype_data_storage(data, index) {}

std::vector<index_t>& cell_state_t::neighbors() { return data_.states.neighbors[index_]; }

std::vector<index_t>& cell_state_t::spring_attachments() { return data_.states.springs[index_]; }

std::vector<index_t>& cell_state_t::attached_cells() { return data_.states.attached_cells[index_]; }

real_t* cell_state_t::orientation() { return data_.states.orientation.data() + index_ * data_.e.mechanics_mesh.dims; }

real_t& cell_state_t::simple_pressure() { return data_.states.simple_pressure[index_]; }

index_t& cell_state_t::number_of_nuclei() { return data_.states.number_of_nuclei[index_]; }

real_t& cell_state_t::damage() { return data_.states.damage[index_]; }

real_t& cell_state_t::total_attack_time() { return data_.states.total_attack_time[index_]; }

cell::cell(agent_id_t id, cell_data& data, index_t index)
	: agent(id, data.agent_data, index), data_(data), phenotype(data_, index), state(data_, index)
{}

real_t* cell::velocity() { return data_.velocities.data() + index_ * data_.e.m.mesh.dims; }

index_t& cell::cell_definition_index() { return data_.cell_definition_indices[index_]; }

uint8_t& cell::is_movable() { return data_.is_movable[index_]; }

void cell::remove() {}
void cell::divide() {}
