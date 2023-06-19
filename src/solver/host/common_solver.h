#pragma once

#include "../../environment.h"

namespace physicell {

class common_solver
{
protected:
	static cell_data& get_cell_data(environment& e);

public:
	static biofvm::index_t get_mesh_index(const biofvm::point_t<biofvm::index_t, 3>& position,
										  const biofvm::cartesian_mesh& mesh);

	static biofvm::point_t<biofvm::index_t, 3> get_mesh_position(const biofvm::real_t* position,
																 const biofvm::cartesian_mesh& mesh);

	template <typename func_t>
	static void for_each_in_mech_neighborhood(const environment& e, const biofvm::point_t<biofvm::index_t, 3>& position,
											  biofvm::index_t i, func_t f)
	{
		// first we process the voxel the agent is in
		// we need to skip all agents that are before the current one because we already processed them
		{
			const auto& cells_in_this_voxel = e.cells_in_mechanics_voxels[get_mesh_index(position, e.mechanics_mesh)];

			auto it = std::find(cells_in_this_voxel.begin(), cells_in_this_voxel.end(), i);

			for (it++; it < cells_in_this_voxel.end(); it++)
				f(*it);
		}

		// now we process the surrounding voxels with increasing indices
		for (biofvm::index_t z = 0; z <= 1; z++)
		{
			if (position[2] + z >= e.mechanics_mesh.grid_shape[2])
				continue;

			for (biofvm::index_t y = 0; y <= 1; y++)
			{
				if (position[1] + y >= e.mechanics_mesh.grid_shape[1])
					continue;

				for (biofvm::index_t x = 0; x <= 1; x++)
				{
					if (position[0] + x >= e.mechanics_mesh.grid_shape[0] || (x == 0 && y == 0 && z == 0))
						continue;

					for (auto& cell_idx : e.cells_in_mechanics_voxels[get_mesh_index(
							 { position[0] + x, position[1] + y, position[2] + z }, e.mechanics_mesh)])
						f(cell_idx);
				}
			}
		}
	}
};

} // namespace physicell
