#pragma once

#include <BioFVM/microenvironment.h>

#include "cell_container.h"
#include "cell_definition.h"
#include "original/user_parameters.h"
#include "runtime_settings.h"

namespace physicell {

struct environment
{
	environment(biofvm::microenvironment& m, biofvm::cartesian_mesh mechanics_mesh);

	biofvm::microenvironment& m;

	cell_container_base& cells();

	template <typename container_t>
	container_t& cast_container()
	{
		return dynamic_cast<container_t&>(*m.agents);
	}

	biofvm::cartesian_mesh mechanics_mesh;

	biofvm::real_t mechanics_time_step;
	biofvm::real_t phenotype_time_step;

	biofvm::real_t current_time;

	runtime_settings settings;
	User_Parameters parameters;

	std::unique_ptr<std::vector<biofvm::index_t>[]> cells_in_mechanics_voxels;

	biofvm::index_t cell_definitions_count;
	std::unique_ptr<cell_definition[]> cell_definitions;
	cell_definition cell_defaults;

	cell_definition* find_cell_definition(std::string name);
	cell_definition* find_cell_definition(biofvm::index_t type);
};

} // namespace physicell
