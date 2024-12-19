#include "src/environment.h"
#include "src/original/modules/settings.h"

namespace physicell {

void make_packed_square(environment& e, cell_definition* cd, biofvm::index_t residency,
						biofvm::point_t<biofvm::real_t, 3> starting_position);

void make_circle(environment& e, cell_definition* cd, biofvm::index_t residency,
				 biofvm::point_t<biofvm::real_t, 3> starting_position);

void setup_potential_parameters(environment& e, User_Parameters& parameters);

} // namespace physicell
