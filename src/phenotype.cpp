#include "phenotype.h"

using namespace biofvm;
using namespace physicell;

phenotype_t::phenotype_t(cell_data& data, const index_t& index)
	: death(data, index),
	  volume(data, index),
	  geometry(data, index),
	  mechanics(data, index),
	  motility(data, index),
	  secretion(data, index),
	  molecular(data, index),
	  cell_interactions(data, index),
	  cell_transformations(data, index),
	  cell_integrity(data, index)
{}

void phenotype_t::copy_to(phenotype_t& dest)
{
	cycle.copy(dest.cycle);
	death.copy(dest.death);
	volume.copy(dest.volume);
	geometry.copy(dest.geometry);
	mechanics.copy(dest.mechanics);
	motility.copy(dest.motility);
	secretion.copy(dest.secretion);
	molecular.copy(dest.molecular);
	cell_interactions.copy(dest.cell_interactions);
	cell_transformations.copy(dest.cell_transformations);
	cell_integrity.copy(dest.cell_integrity);
}

void phenotype_t::copy_to_and_retain_state(phenotype_t& dest)
{
	cycle.copy(dest.cycle);
	death.copy(dest.death);
	
	// copy only subset of volume and no geometry
	dest.volume.target_solid_cytoplasmic() = volume.target_solid_cytoplasmic();
	dest.volume.target_solid_nuclear() = volume.target_solid_nuclear();
	dest.volume.target_fluid_fraction() = volume.target_fluid_fraction();
	dest.volume.target_cytoplasmic_to_nuclear_ratio() = volume.target_cytoplasmic_to_nuclear_ratio();
	dest.volume.relative_rupture_volume() = volume.relative_rupture_volume();

	mechanics.copy(dest.mechanics);
	motility.copy(dest.motility);
	secretion.copy(dest.secretion);
	molecular.copy(dest.molecular);
	cell_interactions.copy(dest.cell_interactions);
	cell_transformations.copy(dest.cell_transformations);
	cell_integrity.copy(dest.cell_integrity);
}