#include "cell_definition.h"

#include "environment.h"

using namespace biofvm;
using namespace physicell;

cell_definition::cell_definition(environment& e, index_t index)
	: data_(e.cell_definitions_data),
	  index(index),
	  type(0),
	  name("unnamed"),
	  is_movable(true),
	  e(e),
	  phenotype(data_, index)
{
	parameters.pReference_live_phenotype = &phenotype;
}

void cell_definition::inherit_from(cell_definition& def)
{
	is_movable = def.is_movable;

	parameters = def.parameters;
	custom_data = def.custom_data;
	functions = def.functions;

	def.phenotype.copy_to(phenotype);

	parameters.pReference_live_phenotype = &phenotype;
}
