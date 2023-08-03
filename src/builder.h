#pragma once

#include <optional>
#include <pugixml.hpp>
#include <string>

#include <BioFVM/microenvironment_builder.h>

#include "environment.h"
#include "original/modules/settings.h"

namespace physicell {

class builder
{
	std::string config_path_;
	std::optional<pugi::xml_document> config_doc_;
	std::optional<pugi::xml_node> config_root_;

	biofvm::microenvironment_builder m_builder_;

	std::optional<PhysiCell_Settings> settings_;
	std::optional<User_Parameters> parameters_;
	std::optional<biofvm::microenvironment> m_;
	std::unique_ptr<environment> e_;

	bool cell_defs_initialized_;
	std::vector<std::string> cell_definition_names_;

	void peek_cell_definitions();

	biofvm::microenvironment& get_microenvironment();

	void construct_single_cell_definition(const pugi::xml_node& node);

	biofvm::index_t find_cell_definition_index(const std::string& name);

	void load_rules();

	void load_signals();

public:
	builder(int argc, char** argv);

	pugi::xml_node& get_config_root();

	// for accessing the settings
	PhysiCell_Settings& get_settings();

	// for accessing the parameters
	User_Parameters& get_parameters();

	// for modifying the microenvironment
	biofvm::microenvironment_builder& get_microenvironment_builder();

	// for modifying default cell definition
	cell_definition& get_default_cell_definition();

	// for modifying cell definitions
	std::vector<cell_definition>& get_cell_definitions();
	cell_definition* find_cell_definition(const std::string& name);

	environment& get_environment();

	std::unique_ptr<environment> build_environment();
};

} // namespace physicell
