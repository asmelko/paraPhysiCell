#include "kelvin_voigt_model.h"

#include <cmath>

#include "../environment.h"

using namespace biofvm;
using namespace physicell;

void kelvin_voigt_model::update_cell_positions(environment& e)
{
	update_cell_forces(e);
	update_motility(e);
	update_basement_membrane_interactions(e);
	update_positions(e);
}
