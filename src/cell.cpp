#include "cell.h"

#include <BioFVM/microenvironment.h>

#include "cell_data.h"
#include "environment.h"
#include "solver/host/containers_solver.h"
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

void cell::set_default(cell_definition& def, index_t def_index)
{
	is_movable() = def.is_movable;

	convert(def, def_index);
}

void cell::convert(cell_definition& def, index_t def_index)
{
	type = def.type;
	type_name = def.name;

	custom_data = def.custom_data;
	parameters = def.parameters;
	functions = def.functions;

	def.phenotype.copy(phenotype);

	cell_definition_index() = def_index;

	assign_orientation();
}

void cell::copy_from(cell& source)
{
	type = source.type;
	type_name = source.type_name;

	custom_data = source.custom_data;
	parameters = source.parameters;
	functions = source.functions;

	phenotype.copy(source.phenotype);

	cell_definition_index() = source.cell_definition_index();

	assign_orientation();

	for (index_t d = 0; d < data_.e.m.mesh.dims; d++)
	{
		velocity()[d] = source.velocity()[d];
	}
}

void cell::divide(cell& new_cell)
{
	// remove all attached
	remove_springs(index_, data_.states.springs.data());
	state.spring_attachments().clear();

	// divide conserved quantitites in custom data in half
	custom_data.divide_conserved_quantities();

	// divide internalized substrates in half
	if (data_.e.m.compute_internalized_substrates)
	{
		phenotype.molecular.divide();
	}

	// clear damage and total attack time
	state.damage() = 0;
	state.total_attack_time() = 0;

	// divide volume in half
	phenotype.volume.divide();
	phenotype.geometry.update();

	// create random vector for position change
	real_t random_walk[3];
	{
		if (data_.e.m.mesh.dims == 1)
			position_helper<1>::random_walk(false, random_walk);
		else if (data_.e.m.mesh.dims == 2)
			position_helper<2>::random_walk(false, random_walk);
		else if (data_.e.m.mesh.dims == 3)
			position_helper<3>::random_walk(false, random_walk);

		real_t tmp = 0;
		for (index_t d = 0; d < data_.e.m.mesh.dims; d++)
		{
			tmp += random_walk[d] * state.orientation()[d];
		}

		for (index_t d = 0; d < data_.e.m.mesh.dims; d++)
		{
			random_walk[d] = (random_walk[d] - tmp * state.orientation()[d] * phenotype.geometry.polarity())
							 * phenotype.geometry.radius();
		}
	}

	new_cell.copy_from(*this);

	for (index_t d = 0; d < data_.e.m.mesh.dims; d++)
	{
		new_cell.position()[d] += random_walk[d];
		position()[d] -= random_walk[d] * 0.5;
	}
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

cell_state_flag& cell::flag() { return data_.flags[index_]; }

void cell::flag_for_removal() { flag() = cell_state_flag::to_remove; }

void cell::flag_for_division() { flag() = cell_state_flag::to_divide; }

const environment& cell::e() const { return data_.e; }
