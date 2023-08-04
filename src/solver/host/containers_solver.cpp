#include "containers_solver.h"

#include <algorithm>

#include "BioFVM/solver.h"
#include "solver_helper.h"

using namespace biofvm;
using namespace physicell;

void containers_solver::update_mechanics_mesh(environment& e)
{
	for (std::size_t i = 0; i < e.mechanics_mesh.voxel_count(); i++)
		e.cells_in_mechanics_voxels[i].clear();

	const auto& data = get_cell_data(e);

	for (index_t i = 0; i < data.agents_count; i++)
	{
		auto mech_pos =
			get_mesh_position(data.agent_data.positions.data() + i * e.mechanics_mesh.dims, e.mechanics_mesh);

		e.cells_in_mechanics_voxels[get_mesh_index(mech_pos, e.mechanics_mesh)].push_back(i);
	}
}

void physicell::remove_attached(index_t to_remove, std::vector<index_t>* __restrict__ attached)
{
	for (auto spring : attached[to_remove])
	{
		auto it = std::find_if(attached[spring].begin(), attached[spring].end(),
							   [to_remove](index_t i) { return i == to_remove; });

		*it = attached[spring].back();
		attached[spring].pop_back();
	}
}

void rename_attached(index_t old_index, index_t new_index, std::vector<index_t>* __restrict__ attached)
{
	for (auto spring : attached[old_index])
	{
		auto it = std::find_if(attached[spring].begin(), attached[spring].end(),
							   [old_index](index_t i) { return i == old_index; });

		*it = new_index;
	}
}

void remove_single(index_t i, const cell_state_flag* __restrict__ flag, const real_t* __restrict__ positions,
				   std::vector<index_t>* __restrict__ springs, std::vector<index_t>* __restrict__ attachced_cells,
				   const cartesian_mesh& mesh, cell_container& container, index_t& counter)
{
	while (true)
	{
		bool out_of_bounds = false;

		for (index_t d = 0; d < mesh.dims; d++)
		{
			if (positions[i * mesh.dims + d] < mesh.bounding_box_mins[d]
				|| positions[i * mesh.dims + d] > mesh.bounding_box_maxs[d])
			{
				out_of_bounds = true;
			}
		}

		if (out_of_bounds == false && flag[i] != cell_state_flag::to_remove)
			return;

		if (flag[i] == cell_state_flag::to_remove)
			counter++;

		remove_attached(i, springs);
		rename_attached(container.data().agents_count - 1, i, springs);

		remove_attached(i, attachced_cells);
		rename_attached(container.data().agents_count - 1, i, attachced_cells);

		container.remove_at(i);
	}
}

void update_cell_container_internal(const cell_state_flag* __restrict__ flag, const real_t* __restrict__ positions,
									std::vector<index_t>* __restrict__ springs,
									std::vector<index_t>* __restrict__ attached_cells, const cartesian_mesh& mesh,
									cell_container& container, index_t& counter)
{
	for (index_t i = 0; i < container.data().agents_count; i++)
	{
		remove_single(i, flag, positions, springs, attached_cells, mesh, container, counter);
	}
}

void containers_solver::update_cell_container_for_mechanics(environment& e)
{
	auto& data = get_cell_data(e);

	update_cell_container_internal(data.flags.data(), data.agent_data.positions.data(), data.states.springs.data(),
								   data.states.attached_cells.data(), e.m.mesh, e.cast_container<cell_container>(),
								   e.deaths_count);
}

void containers_solver::update_cell_container_for_phenotype(environment& e, cell_solver& s)
{
	auto& data = get_cell_data(e);
	auto& c = e.cast_container<cell_container>();

	const auto n = c.agents().size();

	for (std::size_t i = 0; i < n; i++)
	{
		if (c.agents()[i]->flag() == cell_state_flag::to_divide)
		{
			c.agents()[i]->flag() = cell_state_flag::none;
			
			auto cell = c.create();

			c.agents()[i]->divide(*cell);

			e.divisions_count++;
		}
	}

	for (index_t i = 0; i < data.agents_count; i++)
	{
		s.release_internalized_substrates(e.m, i);

		remove_single(i, data.flags.data(), data.agent_data.positions.data(), data.states.springs.data(),
					  data.states.attached_cells.data(), e.m.mesh, c, e.deaths_count);
	}
}
