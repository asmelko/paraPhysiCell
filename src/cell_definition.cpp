#include "cell_definition.h"

using namespace physicell;

Cell_Definition::Cell_Definition(environment& e)
	: data_(e), type(0), name("unnamed"), is_movable(true), e(e), phenotype(data_, 0)
{
	data_.add();
	parameters.pReference_live_phenotype = &phenotype;

	// set up the default functions
	// functions.volume_update_function = NULL; // standard_volume_update_function;
	// functions.update_migration_bias = NULL;

	// functions.update_phenotype = NULL;
	// functions.custom_cell_rule = NULL;

	// functions.update_velocity = NULL; // standard_update_cell_velocity;
	// functions.add_cell_basement_membrane_interactions = NULL;
	// functions.calculate_distance_to_membrane = NULL;

	// functions.set_orientation = NULL;
}
