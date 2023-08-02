#pragma once

#include <BioFVM/microenvironment.h>

#include "cell_container.h"
#include "cell_definition.h"
#include "types.h"

namespace physicell {

constexpr biofvm::point_t<biofvm::index_t, 3> default_mechanics_voxel_shape = { 30, 30, 30 };

struct environment
{
	environment(biofvm::microenvironment& m, biofvm::index_t cell_definitions_count,
				biofvm::point_t<biofvm::index_t, 3> mechanics_voxel_shape = default_mechanics_voxel_shape);

	biofvm::microenvironment& m;

	cell_container_base& cells();

	template <typename container_t>
	container_t& cast_container()
	{
		return dynamic_cast<container_t&>(*m.agents);
	}

	bool virtual_wall_at_domain_edges;
	bool rules_enabled;
	bool automated_spring_adhesion;

	biofvm::cartesian_mesh mechanics_mesh;

	biofvm::real_t mechanics_time_step;
	biofvm::real_t phenotype_time_step;

	biofvm::real_t current_time;

	std::unique_ptr<std::vector<biofvm::index_t>[]> cells_in_mechanics_voxels;

	biofvm::index_t cell_definitions_count;
	std::vector<cell_definition> cell_definitions;
	cell_definition& cell_defaults();

	cell_definition* find_cell_definition(const std::string& name);
	cell_definition* find_cell_definition(biofvm::index_t type);

	void display_info();

	void display_cell_definitions_info();
};

} // namespace physicell
