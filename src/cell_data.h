#pragma once

#include <cstdint>
#include <functional>
#include <vector>

#include <BioFVM/agent_data.h>

namespace physicell {

struct environment;
class cell;

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

	std::vector<biofvm::real_t> cytoplasmic_biomass_change_rate;
	std::vector<biofvm::real_t> nuclear_biomass_change_rate;
	std::vector<biofvm::real_t> fluid_change_rate;
	std::vector<biofvm::real_t> calcification_rate;
	std::vector<biofvm::real_t> target_solid_cytoplasmic;
	std::vector<biofvm::real_t> target_solid_nuclear;
	std::vector<biofvm::real_t> target_fluid_fraction;
	std::vector<biofvm::real_t> target_cytoplasmic_to_nuclear_ratio;
	std::vector<biofvm::real_t> relative_rupture_volume;

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

	void add(biofvm::index_t size, biofvm::index_t cell_definitions_count);
	void remove(biofvm::index_t index, biofvm::index_t size, biofvm::index_t cell_definitions_count);
};

struct motility_data
{
	using direction_update_func = std::function<void(cell&)>;

	std::vector<std::uint8_t> is_motile;
	std::vector<biofvm::real_t> persistence_time;
	std::vector<biofvm::real_t> migration_speed;

	std::vector<biofvm::real_t> migration_bias_direction;
	std::vector<biofvm::real_t> migration_bias;

	std::vector<biofvm::real_t> motility_vector;

	std::vector<std::uint8_t> restrict_to_2d;

	std::vector<biofvm::index_t> chemotaxis_index;
	std::vector<biofvm::index_t> chemotaxis_direction;
	std::vector<biofvm::real_t> chemotactic_sensitivities;

	std::vector<direction_update_func> update_migration_bias_direction;

	void add(biofvm::index_t size, biofvm::index_t dims, biofvm::index_t substrates_count);
	void remove(biofvm::index_t index, biofvm::index_t size, biofvm::index_t dims, biofvm::index_t substrates_count);
};

struct interactions_data
{
	std::vector<biofvm::real_t> dead_phagocytosis_rate;
	std::vector<biofvm::real_t> live_phagocytosis_rates;

	std::vector<biofvm::real_t> damage_rate;
	std::vector<biofvm::real_t> attack_rates;
	std::vector<biofvm::real_t> immunogenicities;

	std::vector<biofvm::real_t> fusion_rates;

	void add(biofvm::index_t size, biofvm::index_t cell_definitions_count);
	void remove(biofvm::index_t index, biofvm::index_t size, biofvm::index_t cell_definitions_count);
};

struct transformations_data
{
	std::vector<biofvm::real_t> transformation_rates;

	void add(biofvm::index_t size, biofvm::index_t cell_definitions_count);
	void remove(biofvm::index_t index, biofvm::index_t size, biofvm::index_t cell_definitions_count);
};

struct death_data
{
	std::vector<std::uint8_t> dead;

	void add(biofvm::index_t size);
	void remove(biofvm::index_t index, biofvm::index_t size);
};

struct cell_state_data
{
	std::vector<std::vector<biofvm::index_t>> neighbors;

	std::vector<std::vector<biofvm::index_t>> springs;
	std::vector<std::vector<biofvm::index_t>> attached_cells;

	std::vector<biofvm::real_t> orientation;
	std::vector<biofvm::real_t> simple_pressure;
	std::vector<biofvm::index_t> number_of_nuclei;

	std::vector<biofvm::real_t> damage;
	std::vector<biofvm::real_t> total_attack_time;

	void add(biofvm::index_t size, biofvm::index_t dims);
	void remove(biofvm::index_t index, biofvm::index_t size, biofvm::index_t dims);
};

enum class cell_state_flag : std::uint8_t
{
	none = 0,
	to_remove = 1,
	to_divide = 2
};

struct cell_data
{
	// BioFVM phenotype data: secretion + total volume + molecular
	biofvm::agent_data agent_data;

	// PhysiCell phenotype data
	volume_data volumes;
	geometry_data geometries;
	mechanics_data mechanics;
	motility_data motilities;
	death_data deaths;

	interactions_data interactions;
	transformations_data transformations;

	cell_state_data states;

	std::vector<biofvm::real_t> previous_velocities;
	std::vector<biofvm::real_t> velocities;
	std::vector<biofvm::index_t> cell_definition_indices;
	std::vector<std::uint8_t> is_movable;

	std::vector<cell_state_flag> flags;

	// references agent_data.agents_count
	biofvm::index_t& agents_count;

	environment& e;

	cell_data(environment& e);
	cell_data(environment& e, biofvm::index_t agents_count);

	void add();
	void remove(biofvm::index_t index);
};

} // namespace physicell
