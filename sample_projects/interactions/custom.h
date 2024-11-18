#pragma once

#include "src/builder.h"
#include "src/original/modules/pathology.h"

using namespace biofvm;
using namespace physicell;

// setup functions to help us along

void create_cell_types(builder& builder);
void setup_tissue(environment& e, User_Parameters& parameters, const pugi::xml_node& config_root);

// set up the BioFVM microenvironment
void setup_microenvironment(microenvironment_builder& m_builder);

// custom pathology coloring function

std::vector<std::string> my_coloring_function(cell*);

std::string my_coloring_function_for_substrate(biofvm::real_t concentration, biofvm::real_t max_conc,
											   biofvm::real_t min_conc, PhysiCell_Settings& settings);
