#include "phenotype.h"

using namespace biofvm;
using namespace physicell;

phenotype_t::phenotype_t(cell_data& data, index_t index)
	: death(data, index),
	  volume(data, index),
	  geometry(data, index),
	  mechanics(data, index),
	  motility(data, index),
	  secretion(data, index),
	  molecular(data, index),
	  interactions(data, index),
	  transformations(data, index)
{}

void phenotype_t::copy(phenotype_t& dest)
{
	cycle.copy(dest.cycle);
	death.copy(dest.death);
	volume.copy(dest.volume);
	geometry.copy(dest.geometry);
	mechanics.copy(dest.mechanics);
	motility.copy(dest.motility);
	secretion.copy(dest.secretion);
	molecular.copy(dest.molecular);
	interactions.copy(dest.interactions);
	transformations.copy(dest.transformations);
}
