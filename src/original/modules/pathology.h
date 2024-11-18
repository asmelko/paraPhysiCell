#pragma once

#include <functional>
#include <string>
#include <vector>

namespace physicell {

struct environment;
class cell;
class cell_definition;
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

using cell_coloring_funct_t = std::function<std::vector<std::string>(cell*)>;
using substrate_coloring_funct_t =
	std::function<std::string(double concentration, double max_conc, double min_conc, PhysiCell_Settings&)>;
using cell_counts_funct_t = std::function<void(char*)>;

void create_plot_legend(std::string filename, cell_coloring_funct_t cell_coloring_function, environment& e);

std::vector<std::string> paint_by_number_cell_coloring(cell* pCell); // done

std::string paint_by_density_percentage(double concentration, double max_conc, double min_conc,
										PhysiCell_Settings& settings); // done

void setup_svg_substrate_colormap(std::vector<std::string>& colormap, PhysiCell_Settings& settings);

void SVG_plot(std::string filename, environment& e, PhysiCell_Settings& settings, double z_slice, double time,
			  cell_coloring_funct_t cell_coloring_function,
			  substrate_coloring_funct_t substrate_coloring_function = paint_by_density_percentage,
			  cell_counts_funct_t cell_counts_function = nullptr);

void standard_agent_SVG(std::ofstream& os, cell* pCell, double z_slice, cell_coloring_funct_t cell_coloring_function,
						double X_lower, double Y_lower);
void standard_agent_legend(std::ofstream& os, cell_definition* cell_definition, double& cursor_x, double& cursor_y,
						   cell_coloring_funct_t cell_coloring_function, double temp_cell_radius, environment& e);

} // namespace physicell
