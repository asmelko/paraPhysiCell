#pragma once

#include "data_storage_phenotype.h"
#include "original/implicit_storage_phenotype.h"

namespace physicell {

struct phenotype_t
{
	Cycle cycle;
	Death death;

	volume_t volume;
	geometry_t geometry;
	mechanics_t mechanics;
	motility_t motility;
	secretion_t secretion;

	molecular_t molecular;

	interactions_t cell_interactions;
	transformations_t cell_transformations;

	phenotype_t(cell_data& data, biofvm::index_t index);

	void copy_to(phenotype_t& dest);
};
} // namespace physicell
