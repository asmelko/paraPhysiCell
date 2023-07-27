#pragma once

#include <pugixml.hpp>
#include <string>

namespace physicell {

struct runtime_settings
{
	bool rules_enabled;
	std::string folder;

	pugi::xml_node config_dom_root;
};

} // namespace physicell
