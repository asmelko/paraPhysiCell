#pragma once

#include <cstdint>
#include <functional>

#include <BioFVM/types.h>

namespace physicell {

class cell;
struct cell_data;

struct phenotype_data_storage
{
protected:
	cell_data& data_;
	biofvm::index_t index_;

public:
	phenotype_data_storage(cell_data& data, biofvm::index_t index);
};

struct volume_t : public phenotype_data_storage
{
	volume_t(cell_data& data, biofvm::index_t index);

	biofvm::real_t& total();

	biofvm::real_t& solid();
	biofvm::real_t& fluid();
	biofvm::real_t& fluid_fraction();

	biofvm::real_t& nuclear();
	biofvm::real_t& nuclear_fluid();
	biofvm::real_t& nuclear_solid();

	biofvm::real_t& cytoplasmic();
	biofvm::real_t& cytoplasmic_fluid();
	biofvm::real_t& cytoplasmic_solid();

	biofvm::real_t& calcified_fraction();

	biofvm::real_t& cytoplasmic_to_nuclear_ratio();

	biofvm::real_t& rupture_volume();

	biofvm::real_t& cytoplasmic_biomass_change_rate();
	biofvm::real_t& nuclear_biomass_change_rate();
	biofvm::real_t& fluid_change_rate();
	biofvm::real_t& calcification_rate();
	biofvm::real_t& target_solid_cytoplasmic();
	biofvm::real_t& target_solid_nuclear();
	biofvm::real_t& target_fluid_fraction();
	biofvm::real_t& target_cytoplasmic_to_nuclear_ratio();
	biofvm::real_t& relative_rupture_volume();

	void divide();

	void multiply_by_factor(biofvm::real_t factor);

	void copy(volume_t& dest);

	void set_defaults();
};

struct geometry_t : public phenotype_data_storage
{
	geometry_t(cell_data& data, biofvm::index_t index);

	biofvm::real_t& radius();
	biofvm::real_t& nuclear_radius();
	biofvm::real_t& surface_area();
	biofvm::real_t& polarity();

	void copy(geometry_t& dest);

	void update();

	void set_defaults();
};

struct mechanics_t : public phenotype_data_storage
{
	mechanics_t(cell_data& data, biofvm::index_t index);

	biofvm::real_t& cell_cell_adhesion_strength();
	biofvm::real_t& cell_BM_adhesion_strength();

	biofvm::real_t& cell_cell_repulsion_strength();
	biofvm::real_t& cell_BM_repulsion_strength();

	biofvm::real_t* cell_adhesion_affinities();

	biofvm::real_t& relative_maximum_adhesion_distance();

	biofvm::index_t& maximum_number_of_attachments();
	biofvm::real_t& attachment_elastic_constant();

	biofvm::real_t& attachment_rate();
	biofvm::real_t& detachment_rate();

	void copy(mechanics_t& dest);

	void set_defaults();
};

struct motility_t : public phenotype_data_storage
{
	motility_t(cell_data& data, biofvm::index_t index);

	std::uint8_t& is_motile();
	biofvm::real_t& persistence_time();
	biofvm::real_t& migration_speed();

	biofvm::real_t* migration_bias_direction();
	biofvm::real_t& migration_bias();

	biofvm::real_t* motility_vector();

	biofvm::index_t& chemotaxis_index();
	biofvm::index_t& chemotaxis_direction();
	biofvm::real_t* chemotactic_sensitivities();

	std::function<void(cell&)>& update_migration_bias_direction();

	void copy(motility_t& dest);

	void set_defaults();
};

struct secretion_t : public phenotype_data_storage
{
	secretion_t(cell_data& data, biofvm::index_t index);

	biofvm::real_t* secretion_rates();
	biofvm::real_t* uptake_rates();
	biofvm::real_t* saturation_densities();
	biofvm::real_t* net_export_rates();

	void copy(secretion_t& dest);

	void set_all_secretion_to_zero();
	void scale_all_uptake_by_factor(biofvm::real_t factor);

	void set_defaults();
};

struct molecular_t : public phenotype_data_storage
{
	molecular_t(cell_data& data, biofvm::index_t index);

	biofvm::real_t* internalized_total_substrates();
	biofvm::real_t* fraction_released_at_death();
	biofvm::real_t* fraction_transferred_when_ingested();

	void divide();

	void copy(molecular_t& dest);

	void set_defaults();
};

struct interactions_t : public phenotype_data_storage
{
	interactions_t(cell_data& data, biofvm::index_t index);

	biofvm::real_t& dead_phagocytosis_rate();
	biofvm::real_t* live_phagocytosis_rates();

	biofvm::real_t& damage_rate();
	biofvm::real_t* attack_rates();
	biofvm::real_t* immunogenicities();

	biofvm::real_t* fusion_rates();

	void copy(interactions_t& dest);

	void set_defaults();
};

struct transformations_t : public phenotype_data_storage
{
	transformations_t(cell_data& data, biofvm::index_t index);

	biofvm::real_t* transformation_rates();

	void copy(transformations_t& dest);

	void set_defaults();
};

} // namespace physicell
