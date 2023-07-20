#include "mechanics_solver.h"

#include <algorithm>

using namespace biofvm;
using namespace physicell;

void mechanics_solver::update_mechanics_mesh(environment& e)
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

void remove_springs(index_t to_remove, std::vector<index_t>* __restrict__ springs)
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

void update_cell_container_internal(index_t n, const std::uint8_t* __restrict__ to_remove,
									const real_t* __restrict__ positions, std::vector<index_t>* __restrict__ springs,
									cell_data& data)
{
	for (index_t i = 0; i < n; i++)
	{
		bool out_of_bounds = false;

		for (index_t d = 0; d < data.e.mechanics_mesh.dims; d++)
		{
			if (positions[i * d + d] < data.e.mechanics_mesh.bounding_box_mins[d]
				|| positions[i * d + d] > data.e.mechanics_mesh.bounding_box_maxs[d])
			{
				out_of_bounds = true;
			}
		}

		if (out_of_bounds == false && to_remove[i] == 0)
			continue;

		remove_springs(i, springs);
		rename_springs(n - 1, i, springs);

		data.remove(i);
	}
}

void mechanics_solver::update_cell_container(environment& e)
{
	auto& data = get_cell_data(e);

	for (index_t i = 0; i < data.agents_count; i++)
	{
		if (data.to_remove[i] == 0)
			continue;

		remove_springs(i, data.states.springs.data());
		rename_springs(data.agents_count - 1, i, data.states.springs.data());

		data.remove(i);
	}
}
