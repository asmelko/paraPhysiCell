#include "cell_definition.h"

using namespace physicell;

Cell_Definition::Cell_Definition(environment& e)
	: data_(e), type(0), name("unnamed"), is_movable(true), e(e), phenotype(data_, 0)
{
	data_.add();
	parameters.pReference_live_phenotype = &phenotype;
}
