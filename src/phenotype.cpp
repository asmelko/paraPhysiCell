#include "phenotype.h"

#include "cell.h"
#include "cell_data.h"
#include "environment.h"

using namespace biofvm;
using namespace physicell;

phenotype_data_storage::phenotype_data_storage(cell_data& data, index_t index) : data_(data), index_(index) {}

volume_t::volume_t(cell_data& data, index_t index) : phenotype_data_storage(data, index) {}

real_t& volume_t::total() { return data_.agent_data.volumes[index_]; }

real_t& volume_t::solid() { return data_.volumes.solid[index_]; }

real_t& volume_t::fluid() { return data_.volumes.fluid[index_]; }

real_t& volume_t::fluid_fraction() { return data_.volumes.fluid_fraction[index_]; }

real_t& volume_t::nuclear() { return data_.volumes.nuclear[index_]; }

real_t& volume_t::nuclear_fluid() { return data_.volumes.nuclear_fluid[index_]; }

real_t& volume_t::nuclear_solid() { return data_.volumes.nuclear_solid[index_]; }

real_t& volume_t::cytoplasmic() { return data_.volumes.cytoplasmic[index_]; }

real_t& volume_t::cytoplasmic_fluid() { return data_.volumes.cytoplasmic_fluid[index_]; }

real_t& volume_t::cytoplasmic_solid() { return data_.volumes.cytoplasmic_solid[index_]; }

real_t& volume_t::calcified_fraction() { return data_.volumes.calcified_fraction[index_]; }

real_t& volume_t::cytoplasmic_to_nuclear_ratio() { return data_.volumes.cytoplasmic_to_nuclear_ratio[index_]; }

real_t& volume_t::rupture_volume() { return data_.volumes.rupture_volume[index_]; }

geometry_t::geometry_t(cell_data& data, index_t index) : phenotype_data_storage(data, index) {}

real_t& geometry_t::radius() { return data_.geometries.radius[index_]; }

real_t& geometry_t::nuclear_radius() { return data_.geometries.nuclear_radius[index_]; }

real_t& geometry_t::surface_area() { return data_.geometries.surface_area[index_]; }

real_t& geometry_t::polarity() { return data_.geometries.polarity[index_]; }

mechanics_t::mechanics_t(cell_data& data, index_t index) : phenotype_data_storage(data, index) {}

real_t& mechanics_t::cell_cell_adhesion_strength() { return data_.mechanics.cell_cell_adhesion_strength[index_]; }

real_t& mechanics_t::cell_BM_adhesion_strength() { return data_.mechanics.cell_BM_adhesion_strength[index_]; }

real_t& mechanics_t::cell_cell_repulsion_strength() { return data_.mechanics.cell_cell_repulsion_strength[index_]; }

real_t& mechanics_t::cell_BM_repulsion_strength() { return data_.mechanics.cell_BM_repulsion_strength[index_]; }

real_t* mechanics_t::cell_adhesion_affinities()
{
	return data_.mechanics.cell_adhesion_affinities.data() + index_ * data_.e.cell_definitions_count;
}

real_t& mechanics_t::relative_maximum_adhesion_distance()
{
	return data_.mechanics.relative_maximum_adhesion_distance[index_];
}

index_t& mechanics_t::maximum_number_of_attachments() { return data_.mechanics.maximum_number_of_attachments[index_]; }

real_t& mechanics_t::attachment_elastic_constant() { return data_.mechanics.attachment_elastic_constant[index_]; }

real_t& mechanics_t::attachment_rate() { return data_.mechanics.attachment_rate[index_]; }

real_t& mechanics_t::detachment_rate() { return data_.mechanics.detachment_rate[index_]; }

phenotype_t::phenotype_t(cell_data& data, index_t index)
	: volume(data, index), geometry(data, index), mechanics(data, index)
{}
