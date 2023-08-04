#pragma once

#include <functional>
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
			  std::function<std::vector<std::string>(cell*)> cell_coloring_function,
			  std::function<std::vector<std::string>(double, double, double)> substrate_coloring_function = nullptr);

void create_plot_legend(std::string filename, std::function<std::vector<std::string>(cell*)> cell_coloring_function,
						environment& e);

std::vector<std::string> paint_by_number_cell_coloring(cell* pCell); // done

} // namespace physicell
