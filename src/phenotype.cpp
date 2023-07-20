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
