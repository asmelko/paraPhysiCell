#include "environment.h"

using namespace biofvm;
using namespace physicell;

environment::environment(microenvironment& m, cartesian_mesh mechanics_mesh) : m(m), mechanics_mesh(std::move(mechanics_mesh)) 
{
    cells_in_mechanics_voxels = std::make_unique<std::vector<index_t>[]>(mechanics_mesh.voxel_count());
}

cell_container_base& environment::cells() { return cast_container<cell_container_base&>(); }
