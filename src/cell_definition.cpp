#include "cell_definition.h"

using namespace physicell;

cell_definition::cell_definition(environment& e)
	: data_(e), type(0), name("unnamed"), is_movable(true), e(e), phenotype(data_, 0)
{
	data_.add();
	parameters.pReference_live_phenotype = &phenotype;
}
