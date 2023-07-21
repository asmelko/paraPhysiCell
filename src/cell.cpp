#include "cell.h"

#include <BioFVM/microenvironment.h>

#include "cell_data.h"
#include "environment.h"
#include "solver/host/solver_helper.h"

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

void cell::set_default(cell_definition& def)
{
	is_movable() = def.is_movable;

	convert(def);
}

void cell::convert(cell_definition& def)
{
	type = def.type;
	type_name = def.name;

	custom_data = def.custom_data;
	parameters = def.parameters;
	functions = def.functions;

	def.phenotype.copy(phenotype);

	assign_orientation();
}

void cell::assign_orientation()
{
	if (functions.set_orientation)
	{
		functions.set_orientation(*this, 0);
	}
	else
	{
		if (data_.e.mechanics_mesh.dims == 1)
		{
			position_helper<1>::random_walk(false, state.orientation());
		}
		else if (data_.e.mechanics_mesh.dims == 2)
		{
			position_helper<2>::random_walk(false, state.orientation());
		}
		else if (data_.e.mechanics_mesh.dims == 3)
		{
			position_helper<3>::random_walk(false, state.orientation());
		}
	}
}

void cell::remove() {}
void cell::divide() {}
