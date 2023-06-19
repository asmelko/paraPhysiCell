#include "common_solver.h"

#include <noarr/structures/extra/shortcuts.hpp>

#include "noarr/structures/extra/funcs.hpp"
#include "types.h"

using namespace biofvm;
using namespace physicell;

cell_data& common_solver::get_cell_data(environment& e) { return e.cells().get_cell_data(); }

index_t common_solver::get_mesh_index(const point_t<index_t, 3>& position, const cartesian_mesh& mesh)
{
	auto mesh_l = noarr::scalar<uint8_t>()
				  ^ noarr::sized_vectors<'x', 'y', 'z'>(mesh.grid_shape[0], mesh.grid_shape[1], mesh.grid_shape[2]);

	return mesh_l | noarr::offset<'x', 'y', 'z'>(position[0], position[1], position[2]);
}

point_t<index_t, 3> common_solver::get_mesh_position(const real_t* position, const biofvm::cartesian_mesh& mesh)
{
	if (mesh.dims == 1)
	{
		point_t<float, 3> pos { position[0], 0, 0 };
		auto voxel_pos = mesh.voxel_position(pos);
		return { voxel_pos[0], 0, 0 };
	}
	else if (mesh.dims == 2)
	{
		point_t<float, 3> pos { position[0], position[1], 0 };
		auto voxel_pos = mesh.voxel_position(pos);
		return { voxel_pos[0], voxel_pos[1], 0 };
	}
	else
	{
		point_t<float, 3> pos { position[0], position[1], position[2] };
		auto voxel_pos = mesh.voxel_position(pos);
		return voxel_pos;
	}
}
