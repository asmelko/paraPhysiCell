#include "mechanics_solver.h"

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
