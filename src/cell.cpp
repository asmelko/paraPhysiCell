#include "cell.h"

#include <BioFVM/microenvironment.h>

#include "cell_data.h"
#include "environment.h"
#include "solver/host/containers_solver.h"
#include "solver/host/solver_helper.h"

using namespace biofvm;
using namespace physicell;

cell_state_t::cell_state_t(cell_data& data, const index_t& index) : phenotype_data_storage(data, index) {}

std::vector<index_t>& cell_state_t::neighbors() { return data_.states.neighbors[index_]; }

std::vector<index_t>& cell_state_t::spring_attachments() { return data_.states.springs[index_]; }

std::vector<index_t>& cell_state_t::attached_cells() { return data_.states.attached_cells[index_]; }

real_t* cell_state_t::orientation() { return data_.states.orientation.data() + index_ * data_.e.mechanics_mesh.dims; }

real_t& cell_state_t::simple_pressure() { return data_.states.simple_pressure[index_]; }

index_t& cell_state_t::number_of_nuclei() { return data_.states.number_of_nuclei[index_]; }

real_t& cell_state_t::damage() { return data_.states.damage[index_]; }

real_t& cell_state_t::total_attack_time() { return data_.states.total_attack_time[index_]; }

void cell_state_t::set_defaults()
{
	neighbors().clear();
	spring_attachments().clear();
	attached_cells().clear();

	std::fill(orientation(), orientation() + data_.e.mechanics_mesh.dims, 0);
	simple_pressure() = 0;

	number_of_nuclei() = 1;
	damage() = 0;
	total_attack_time() = 0;
}

cell::cell(agent_id_t id, cell_data& data, index_t index)
	: agent(id, data.agent_data, index), data_(data), phenotype(data_, index_), state(data_, index_)
{
	std::fill(data_.previous_velocities.data() + index_ * data_.e.mechanics_mesh.dims,
			  data_.previous_velocities.data() + (index_ + 1) * data_.e.mechanics_mesh.dims, 0.0);

	std::fill(data_.velocities.data() + index_ * data_.e.mechanics_mesh.dims,
			  data_.velocities.data() + (index_ + 1) * data_.e.mechanics_mesh.dims, 0.0);

	flag() = cell_state_flag::none;

	is_movable() = true;

	state.set_defaults();
}

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

	def.phenotype.copy_to(phenotype);

	cell_definition_index() = def.index;

	assign_orientation();
}

void cell::copy_from(cell& source)
{
	type = source.type;
	type_name = source.type_name;

	custom_data = source.custom_data;
	parameters = source.parameters;
	functions = source.functions;

	source.phenotype.copy_to(phenotype);

	cell_definition_index() = source.cell_definition_index();

	assign_orientation();
}

void cell::divide(cell& new_cell)
{
	// remove all attached
	remove_attached(index_, data_.states.springs.data());
	remove_attached(index_, data_.states.attached_cells.data());
	state.spring_attachments().clear();
	state.attached_cells().clear();

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

	// divide volume in half
	phenotype.volume.divide();
	phenotype.geometry.update();

	new_cell.copy_from(*this);

	for (index_t d = 0; d < data_.e.m.mesh.dims; d++)
	{
		new_cell.position()[d] = position()[d] + random_walk[d];
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

environment& cell::e() { return data_.e; }

cell_container& cell::container() { return data_.e.cast_container<cell_container>(); }

void cell::attach_cell(cell& cell)
{
	auto it = std::find_if(state.attached_cells().begin(), state.attached_cells().end(),
						   [to_find = cell.index()](const index_t index) { return index == to_find; });
	if (it == state.attached_cells().end())
	{
		state.attached_cells().push_back(cell.index());
	}
}

void cell::detach_cell(cell& cell)
{
	auto it = std::find_if(state.attached_cells().begin(), state.attached_cells().end(),
						   [to_find = cell.index()](const index_t index) { return index == to_find; });
	if (it != state.attached_cells().end())
	{
		*it = state.attached_cells().back();
		state.attached_cells().pop_back();
	}
}

void cell::attach_cells(cell& lhs, cell& rhs)
{
	lhs.attach_cell(rhs);
	rhs.attach_cell(lhs);
}

void cell::detach_cells(cell& lhs, cell& rhs)
{
	lhs.detach_cell(rhs);
	rhs.detach_cell(lhs);
}

void cell::remove_all_attached_cells()
{
	remove_attached(index_, data_.states.attached_cells.data());
	state.attached_cells().clear();
}

void cell::assign_position(const biofvm::point_t<biofvm::real_t, 3>& new_position)
{
	for (index_t d = 0; d < data_.e.m.mesh.dims; d++)
	{
		position()[d] = new_position[d];
	}
}

void cell::set_total_volume(real_t volume)
{
	// If the new volume is significantly different than the
	// current total volume, adjust all the sub-volumes
	// proportionally.

	// if( fabs( phenotype.volume.total - volume ) < 1e-16 )
	if (std::abs(phenotype.volume.total() - volume) > 1e-16)
	{
		real_t ratio = volume / (phenotype.volume.total() + 1e-16);
		phenotype.volume.multiply_by_factor(ratio);
	}

	phenotype.geometry.update();
}

std::vector<index_t>& cell::cells_in_my_mechanics_voxel()
{
	auto pos = common_solver::get_mesh_position(position(), e().mechanics_mesh);
	auto idx = common_solver::get_mesh_index(pos, e().mechanics_mesh);

	return e().cells_in_mechanics_voxels[idx];
}
