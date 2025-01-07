#pragma once

#include "morse_position_model.h"

namespace physicell {

class morse_kv_membrane_position_model : public morse_position_model
{
protected:
	static void update_cell_forces(environment& e);

public:
	virtual void update_cell_positions(environment& e) override;

	virtual void update_cell_neighbors(environment& e) override;
};

} // namespace physicell
