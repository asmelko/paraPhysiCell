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

void physicell::remove_springs(index_t to_remove, std::vector<index_t>* __restrict__ springs)
{
	for (auto spring : springs[to_remove])
	{
		auto it = std::find_if(springs[spring].begin(), springs[spring].end(),
							   [to_remove](index_t i) { return i == to_remove; });

		*it = springs[spring].back();
		springs[spring].pop_back();
	}
}

void rename_springs(index_t old_index, index_t new_index, std::vector<index_t>* __restrict__ springs)
{
	for (auto spring : springs[old_index])
	{
		auto it = std::find_if(springs[spring].begin(), springs[spring].end(),
							   [old_index](index_t i) { return i == old_index; });

		*it = new_index;
	}
}

void remove_single(index_t i, const cell_state_flag* __restrict__ flag, const real_t* __restrict__ positions,
				   std::vector<index_t>* __restrict__ springs, const cartesian_mesh& mesh, cell_container& container)
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

		if (out_of_bounds == false && flag[i] == cell_state_flag::none)
			return;

		remove_springs(i, springs);
		rename_springs(container.data().agents_count - 1, i, springs);

		container.remove_at(i);
	}
}

void update_cell_container_internal(const cell_state_flag* __restrict__ flag, const real_t* __restrict__ positions,
									std::vector<index_t>* __restrict__ springs, const cartesian_mesh& mesh,
									cell_container& container)
{
	for (index_t i = 0; i < container.data().agents_count; i++)
	{
		remove_single(i, flag, positions, springs, mesh, container);
	}
}

void containers_solver::update_cell_container_for_mechanics(environment& e)
{
	auto& data = get_cell_data(e);

	update_cell_container_internal(data.flags.data(), data.agent_data.positions.data(), data.states.springs.data(),
								   e.m.mesh, e.cast_container<cell_container>());
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
			auto cell = c.create();

			c.agents()[i]->divide(*cell);
		}
	}

	for (index_t i = 0; i < data.agents_count; i++)
	{
		s.release_internalized_substrates(e.m, i);

		remove_single(i, data.flags.data(), data.agent_data.positions.data(), data.states.springs.data(), e.m.mesh, c);
	}
}
