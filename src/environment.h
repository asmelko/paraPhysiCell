#pragma once

#include <BioFVM/microenvironment.h>

#include "cell_container.h"
#include "cell_definition.h"
#include "models/position_model.h"

namespace physicell {

constexpr biofvm::point_t<biofvm::index_t, 3> default_mechanics_voxel_shape = { 30, 30, 30 };

struct environment
{
	environment(biofvm::microenvironment& m, biofvm::index_t cell_definitions_count,
				biofvm::point_t<biofvm::index_t, 3> mechanics_voxel_shape = default_mechanics_voxel_shape);

	// TODO: fix so it does not have to be always unique ptr
	environment(environment&&) = delete;

	biofvm::microenvironment& m;

	cell_container_base& container_base();

	template <typename container_t = cell_container>
	container_t& get_container()
	{
		return dynamic_cast<container_t&>(*m.agents);
	}

	bool virtual_wall_at_domain_edges;
	bool rules_enabled;
	bool automated_spring_adhesion;

	biofvm::index_t divisions_count;
	biofvm::index_t deaths_count;

	biofvm::cartesian_mesh mechanics_mesh;

	biofvm::real_t mechanics_time_step;
	biofvm::real_t phenotype_time_step;

	double current_time;

	std::unique_ptr<std::vector<biofvm::index_t>[]> cells_in_mechanics_voxels;

	biofvm::index_t cell_definitions_count;
	std::vector<std::unique_ptr<cell_definition>> cell_definitions;
	cell_data cell_definitions_data;

	std::unique_ptr<position_model> position;

	// morse 
	std::vector<biofvm::real_t> inter_scaling_factors, inter_equilibrium_distances, inter_stiffnesses;

	cell_definition& create_cell_definition();

	cell_definition& cell_defaults();

	cell_definition* find_cell_definition(const std::string& name);
	cell_definition* find_cell_definition(biofvm::index_t type);

	void display_info();

	void display_cell_definitions_info();
};

} // namespace physicell
