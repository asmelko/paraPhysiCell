#pragma once

#include <pugixml.hpp>
#include <string>

#include "../../environment.h"

namespace physicell {

bool load_cells_from_pugixml(const pugi::xml_node& root, environment& e);

} // namespace physicell
