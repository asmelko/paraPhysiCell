#include "morse_position_model.h"

using namespace physicell;

void morse_position_model::update_cell_positions(environment& e)
{
	update_cell_forces(e);
	update_motility(e);
	update_basement_membrane_interactions(e);
	update_positions(e);
}
