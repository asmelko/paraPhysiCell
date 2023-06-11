#pragma once

#include <BioFVM/microenvironment.h>

#include "cell_container.h"

namespace physicell {

struct environment
{
public:
	environment(biofvm::microenvironment& m);

	biofvm::microenvironment& m;

	cell_container_base& cells();

	template <typename container_t>
	container_t& cast_container()
	{
		return dynamic_cast<container_t&>(*m.agents);
	}

	// biofvm::cartesian_mesh mechanics_mesh;

	// std::unique_ptr<std::vector<biofvm::index_t>[]> cells_in_mechanics_voxels;

	// TODO: cell definitions
	biofvm::index_t cell_definitions_count;
};

} // namespace physicell
