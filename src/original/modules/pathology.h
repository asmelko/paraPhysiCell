#pragma once

#include <string>
#include <vector>

namespace physicell {

struct environment;
class cell;
class PhysiCell_Settings;

struct PhysiCell_SVG_options_struct
{
	bool plot_nuclei = true;

	std::string simulation_time_units = "min";
	std::string mu = "&#956;";
	std::string simulation_space_units = "&#956;m";

	std::string label_time_units = "days";

	double font_size = 200;
	std::string font_color = "black";
	std::string font = "Arial";

	double length_bar = 100;
};

extern PhysiCell_SVG_options_struct PhysiCell_SVG_options;

void SVG_plot(std::string filename, environment& e, const PhysiCell_Settings& settings, double z_slice, double time,
			  std::vector<std::string> (*cell_coloring_function)(cell*),
			  std::vector<std::string> (*substrate_coloring_function)(double, double, double) = nullptr);

void create_plot_legend(std::string filename, std::vector<std::string> (*cell_coloring_function)(cell*),
						environment& e);

std::vector<std::string> paint_by_number_cell_coloring(cell* pCell); // done

} // namespace physicell