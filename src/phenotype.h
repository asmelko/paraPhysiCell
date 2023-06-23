#pragma once

#include <cstdint>

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

struct phenotype_implicit_storage
{};

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
};

struct geometry_t : public phenotype_data_storage
{
	geometry_t(cell_data& data, biofvm::index_t index);

	biofvm::real_t& radius();
	biofvm::real_t& nuclear_radius();
	biofvm::real_t& surface_area();
	biofvm::real_t& polarity();
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
};

struct molecular_t : public phenotype_data_storage
{
	molecular_t(cell_data& data, biofvm::index_t index);

	biofvm::real_t* internalized_total_substrates();
	biofvm::real_t* fraction_released_at_death();
	biofvm::real_t* fraction_transferred_when_ingested();
};

struct phenotype_t
{
	volume_t volume;
	geometry_t geometry;
	mechanics_t mechanics;
	motility_t motility;
	molecular_t molecular;

	phenotype_t(cell_data& data, biofvm::index_t index);
};

} // namespace physicell
