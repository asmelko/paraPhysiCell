#include "cell_definition.h"

using namespace physicell;

cell_definition::cell_definition(environment& e)
	: data_(e, 1), type(0), name("unnamed"), is_movable(true), e(e), phenotype(data_, 0)
{
	parameters.pReference_live_phenotype = &phenotype;
}

cell_definition cell_definition::create_copy()
{
	cell_definition copy(e);

	copy.type = type;
	copy.name = name;
	copy.is_movable = is_movable;

	copy.parameters = parameters;
	copy.custom_data = custom_data;
	copy.functions = functions;

	copy.phenotype.copy(phenotype);

	copy.parameters.pReference_live_phenotype = &copy.phenotype;

	return copy;
}

void cell_definition::inherit_from(cell_definition& def)
{
	is_movable = def.is_movable;

	parameters = def.parameters;
	custom_data = def.custom_data;
	functions = def.functions;

	phenotype.copy(def.phenotype);

	parameters.pReference_live_phenotype = &phenotype;
}
