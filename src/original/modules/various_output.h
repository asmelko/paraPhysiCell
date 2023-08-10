#pragma once

#include <iostream>

namespace physicell {

struct environment;
class PhysiCell_Settings;

void display_simulation_status(std::ostream& os, environment& e, PhysiCell_Settings& settings);

void log_output(double t, int output_index, environment& e, PhysiCell_Settings& settings, std::ofstream& report_file);

}; // namespace physicell
