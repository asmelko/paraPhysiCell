#include "cell_definition.h"

using namespace physicell;

cell_definition::cell_definition(environment& e)
	: data_(e, 1), type(0), name("unnamed"), is_movable(true), e(e), phenotype(data_, 0)
{
	parameters.pReference_live_phenotype = &phenotype;
}
