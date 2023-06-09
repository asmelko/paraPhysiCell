#pragma once

#include <vector>

#include <BioFVM/agent_data.h>

namespace biofvm {
struct microenvironment;
}

namespace physicell {

struct volume_data
{
	std::vector<biofvm::real_t> solid;
	std::vector<biofvm::real_t> fluid;
	std::vector<biofvm::real_t> fluid_fraction;

	std::vector<biofvm::real_t> nuclear;
	std::vector<biofvm::real_t> nuclear_fluid;
	std::vector<biofvm::real_t> nuclear_solid;

	std::vector<biofvm::real_t> cytoplasmic;
	std::vector<biofvm::real_t> cytoplasmic_fluid;
	std::vector<biofvm::real_t> cytoplasmic_solid;

	std::vector<biofvm::real_t> calcified_fraction;

	std::vector<biofvm::real_t> cytoplasmic_to_nuclear_ratio;

	std::vector<biofvm::real_t> rupture_volume;

	void add(biofvm::index_t size);
	void remove(biofvm::index_t index, biofvm::index_t size);
};

struct geometry_data
{
	std::vector<biofvm::real_t> radius;
	std::vector<biofvm::real_t> nuclear_radius;
	std::vector<biofvm::real_t> surface_area;
	std::vector<biofvm::real_t> polarity;

	void add(biofvm::index_t size);
	void remove(biofvm::index_t index, biofvm::index_t size);
};

struct mechanics_data
{
	biofvm::index_t cell_definitions_count;

	std::vector<biofvm::real_t> cell_cell_adhesion_strength;
	std::vector<biofvm::real_t> cell_BM_adhesion_strength;

	std::vector<biofvm::real_t> cell_cell_repulsion_strength;
	std::vector<biofvm::real_t> cell_BM_repulsion_strength;

	std::vector<biofvm::real_t> cell_adhesion_affinities;

	std::vector<biofvm::real_t> relative_maximum_adhesion_distance;

	std::vector<biofvm::index_t> maximum_number_of_attachments;
	std::vector<biofvm::real_t> attachment_elastic_constant;

	std::vector<biofvm::real_t> attachment_rate;
	std::vector<biofvm::real_t> detachment_rate;

	mechanics_data(biofvm::index_t cell_definitions_count);

	void add(biofvm::index_t size);
	void remove(biofvm::index_t index, biofvm::index_t size);
};

struct cell_data
{
	// BioFVM phenotype data: secretion + total volume + molecular
	biofvm::agent_data agent_data;

	// PhysiCell phenotype data
	volume_data volume;
	geometry_data geometry;
	mechanics_data mechanics;

	std::vector<biofvm::real_t> velocities;
	std::vector<biofvm::index_t> cell_definition_index;
	std::vector<biofvm::real_t> simple_pressure;


	// references agent_data.agents_count
	biofvm::index_t& agents_count;

	biofvm::microenvironment& m;

	cell_data(biofvm::microenvironment& m, biofvm::index_t cell_definitions_count);

	void add();
	void remove(biofvm::index_t index);
};

} // namespace physicell
