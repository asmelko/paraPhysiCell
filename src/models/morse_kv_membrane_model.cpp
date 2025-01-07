#include "morse_kv_membrane_model.h"

using namespace physicell;

void morse_kv_membrane_position_model::update_cell_positions(environment& e)
{
	morse_position_model::update_cell_forces(e);
	update_cell_forces(e);
	update_motility(e);
	update_basement_membrane_interactions(e);
	update_positions(e);
}
