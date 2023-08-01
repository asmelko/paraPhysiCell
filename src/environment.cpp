#include "environment.h"

using namespace biofvm;
using namespace physicell;

environment::environment(microenvironment& m, index_t cell_definitions_count,
						 biofvm::point_t<biofvm::index_t, 3> mechanics_voxel_shape)
	: m(m),
	  virtual_wall_at_domain_edges(false),
	  rules_enabled(false),
	  automated_spring_adhesion(true),
	  mechanics_mesh(m.mesh.dims, m.mesh.bounding_box_mins, m.mesh.bounding_box_maxs, mechanics_voxel_shape),
	  mechanics_time_step(0.1),
	  phenotype_time_step(6),
	  current_time(0),
	  cell_definitions_count(cell_definitions_count)
{
	cells_in_mechanics_voxels = std::make_unique<std::vector<index_t>[]>(mechanics_mesh.voxel_count());
}

cell_container_base& environment::cells() { return cast_container<cell_container_base&>(); }

cell_definition& environment::cell_defaults() { return cell_definitions[0]; }
