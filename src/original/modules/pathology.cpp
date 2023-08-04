#include "pathology.h"

#include <cmath>
#include <cstring>

#include "../../cell.h"
#include "../../environment.h"
#include "../constants.h"
#include "settings.h"
#include "svg.h"
#include "timer.h"

using namespace biofvm;

namespace physicell {

PhysiCell_SVG_options_struct PhysiCell_SVG_options;

std::string formatted_minutes_to_DDHHMM(double minutes)
{
	static std::string output;
	output.resize(1024);

	int nMinutes = std::rint(minutes); // round( minutes );
	// int nDays = (int) floor( (minutes+1e-6) / 1440.0 ); // minutes / 1440
	int nDays = nMinutes / 1440;
	nMinutes -= nDays * 1440;

	// int nHours = (int) floor( (nMinutes+1e-6) / 60.0 ); // nMinutes / 60;
	int nHours = nMinutes / 60;
	double dMinutes = minutes - 60 * (nDays * 24 + nHours);
	if (dMinutes < 0)
	{
		dMinutes = 0.0;
	}
	sprintf((char*)output.c_str(), "%d days, %d hours, and %2.2f minutes", nDays, nHours, dMinutes);

	return output;
}

void SVG_plot(std::string filename, environment& e, const PhysiCell_Settings& settings, double z_slice, double time,
			  cell_coloring_funct_t cell_coloring_function,
			  std::function<std::vector<std::string>(double, double, double)> substrate_coloring_function)
{
	double X_lower = e.m.mesh.bounding_box_mins[0];
	double X_upper = e.m.mesh.bounding_box_maxs[0];

	double Y_lower = e.m.mesh.bounding_box_mins[1];
	double Y_upper = e.m.mesh.bounding_box_maxs[1];

	double plot_width = X_upper - X_lower;
	double plot_height = Y_upper - Y_lower;

	double font_size = 0.025 * plot_height; // PhysiCell_SVG_options.font_size;
	double top_margin = font_size * (.2 + 1 + .2 + .9 + .5);

	// open the file, write a basic "header"
	std::ofstream os(filename, std::ios::out);
	if (os.fail())
	{
		std::cout << std::endl << "Error: Failed to open " << filename << " for SVG writing." << std::endl << std::endl;

		std::cout << std::endl
				  << "Error: We're not writing data like we expect. " << std::endl
				  << "Check to make sure your save directory exists. " << std::endl
				  << std::endl
				  << "I'm going to exit with a crash code of -1 now until " << std::endl
				  << "you fix your directory. Sorry!" << std::endl
				  << std::endl;
		exit(-1);
	}

	if (settings.enable_substrate_plot == true && substrate_coloring_function != NULL)
	{
		double legend_padding = 200.0; // I have to add a margin on the left to visualize the bar plot and the values

		Write_SVG_start(os, plot_width + legend_padding, plot_height + top_margin);

		// draw the background
		Write_SVG_rect(os, 0, 0, plot_width + legend_padding, plot_height + top_margin, 0.002 * plot_height, "white",
					   "white");
	}
	else
	{
		Write_SVG_start(os, plot_width, plot_height + top_margin);

		// draw the background
		Write_SVG_rect(os, 0, 0, plot_width, plot_height + top_margin, 0.002 * plot_height, "white", "white");
	}
	// write the simulation time to the top of the plot

	char* szString;
	szString = new char[1024];

	int total_cell_count = e.m.agents->agents_count();

	double temp_time = time;

	std::string time_label = formatted_minutes_to_DDHHMM(temp_time);

	sprintf(szString, "Current time: %s, z = %3.2f %s", time_label.c_str(), z_slice,
			PhysiCell_SVG_options.simulation_space_units.c_str());
	Write_SVG_text(os, szString, font_size * 0.5, font_size * (.2 + 1), font_size,
				   PhysiCell_SVG_options.font_color.c_str(), PhysiCell_SVG_options.font.c_str());
	sprintf(szString, "%u agents", total_cell_count);
	Write_SVG_text(os, szString, font_size * 0.5, font_size * (.2 + 1 + .2 + .9), 0.95 * font_size,
				   PhysiCell_SVG_options.font_color.c_str(), PhysiCell_SVG_options.font.c_str());

	delete[] szString;


	// add an outer "g" for coordinate transforms

	os << " <g id=\"tissue\" " << std::endl
	   << "    transform=\"translate(0," << plot_height + top_margin << ") scale(1,-1)\">" << std::endl;

	// prepare to do mesh-based plot (later)

	double dx_stroma = e.m.mesh.voxel_shape[0];
	double dy_stroma = e.m.mesh.voxel_shape[1];

	os << "  <g id=\"ECM\">" << std::endl;

	// int ratio = 1;
	// double voxel_size = dx_stroma / (double)ratio;

	// double half_voxel_size = voxel_size / 2.0;
	// double normalizer = 78.539816339744831 / (voxel_size * voxel_size * voxel_size);

	// color in the background ECM
	if (settings.enable_substrate_plot == true && substrate_coloring_function != NULL)
	{
		double dz_stroma = e.m.mesh.voxel_shape[2];
		double max_conc;
		double min_conc;

		std::string sub = settings.substrate_to_monitor;
		int sub_index = e.m.find_substrate_index(sub); // check the substrate does actually exist
		if (sub_index == -1)
		{
			std::cout << "ERROR SAMPLING THE SUBSTRATE: COULD NOT FIND THE SUBSTRATE " << sub
					  << std::endl; // if not print error message
		}
		else
		{
			if (settings.limits_substrate_plot)
			{
				max_conc = settings.max_concentration;
				min_conc = settings.min_concentration;
			}
			else
			{
				max_conc = e.m.substrate_densities[5 * e.m.substrates_count + sub_index];
				min_conc = e.m.substrate_densities[5 * e.m.substrates_count
												   + sub_index]; // so here I am sampling the concentration to set a min
																 // and a mx
				// look for the max and min concentration among all the substrates
				for (std::size_t n = 0; n < e.m.mesh.voxel_count(); n++)
				{
					double concentration = e.m.substrate_densities[n * e.m.substrates_count + sub_index];
					if (concentration > max_conc)
						max_conc = concentration;
					if (concentration < min_conc)
						min_conc = concentration;
				}
			};

			// check that max conc is not zero otherwise it is a big problem!
			if (max_conc == 0)
			{
				max_conc = 1.0;
			};

			for (index_t z = 0; z < e.m.mesh.grid_shape[2]; z++)
				for (index_t y = 0; y < e.m.mesh.grid_shape[1]; y++)
					for (index_t x = 0; x < e.m.mesh.grid_shape[0]; x++)
					{
						auto voxel_center = e.m.mesh.voxel_center({ x, y, z });
						int z_center = voxel_center[2];
						double z_displ = z_center - dz_stroma / 2;

						double z_compare = z_displ;

						if (e.m.mesh.dims == 2)
						{
							z_compare = z_center;
						};

						if (z_slice == z_compare)
						{ // this is to make sure the substrate is sampled in the voxel visualized (so basically the
						  // slice)
							int x_center = voxel_center[0];
							int y_center = voxel_center[1];

							double x_displ = x_center - dx_stroma / 2;
							double y_displ = (y_center - dy_stroma) + dy_stroma / 2;

							double concentration = e.m.substrate_density_value({ x, y, z }, sub_index);

							std::vector<std::string> output =
								substrate_coloring_function(concentration, max_conc, min_conc);

							Write_SVG_rect(os, x_displ - X_lower, y_displ - Y_lower, dx_stroma, dy_stroma, 0, "none",
										   output[0]);
						}
					}

			// add legend for the substrate

			os << " <g id=\"legend\" " << std::endl
			   << "    transform=\"translate(0," << plot_height + 25 << ") scale(1,-1)\">"
			   << std::endl; // for some misterious reasons, the tissue part in the SVG is rotated so I have to
							 // re-rotate to draw
							 //  the legend, otherwise it will be printed upside down
			// int padding = 0;

			double conc_interval =
				(max_conc - min_conc) / 13; // setting the interval for the values in the legend. I will divide the
											// legend in 13 parts (as in the jupyter notebook)

			for (int i = 0; i <= 12; i++)
			{ // creating 13 rectangoles for the bar, each one with a different shade of color.

				double concentration_sample =
					min_conc + (conc_interval * i); // the color depends on the concentration, starting from the min
													// concentration to the max (which was sampled before)

				std::vector<std::string> output = substrate_coloring_function(concentration_sample, max_conc, min_conc);

				// padding = 25 * i;

				double upper_left_x = plot_width + 25.0;
				double upper_left_y = ((plot_height - 25) / 13.0) * i; // here I set the position of each rectangole

				Write_SVG_rect(os, upper_left_x, plot_height - upper_left_y - 60, 25.0, ((plot_height - 25.0) / 13.0),
							   0, "none", output[0]); // drawing each piece of the barplot

				if (i % 2 == 0)
				{ // of course I am not printing each value of the barplot, otherwise is too crowded, so just one each 2

					char* szString;
					szString = new char[1024];

					sprintf(szString, "- %e", concentration_sample);

					Write_SVG_text(os, szString, upper_left_x + 24, plot_height - upper_left_y + 5.31, font_size,
								   PhysiCell_SVG_options.font_color.c_str(),
								   PhysiCell_SVG_options.font
									   .c_str()); // misterious values set with a trial and error approach due to OCD.
												  // But now the legend is coherent at pixel level

					delete[] szString;
				}
			}

			Write_SVG_rect(os, 25.0 + plot_width, 25.0, 25.0, plot_height - 25.0, 0.002 * plot_height, "black",
						   "none"); // nice black contour around the legend

			os << "  </g>" << std::endl; // no more rotation, restoring the tissue object in the SVG
		}
	}
	/*
	 if( ECM.TellRows() > 0 )
	 {
	  // find the k corresponding to z_slice



	  Vector position;
	  *position(2) = z_slice;


	  // 25*pi* 5 microns^2 * length (in source) / voxelsize^3

	  for( int j=0; j < ratio*ECM.TellCols() ; j++ )
	  {
	   // *position(1) = *Y_environment(j);
	   *position(1) = *Y_environment(0) - dy_stroma/2.0 + j*voxel_size + half_voxel_size;

	   for( int i=0; i < ratio*ECM.TellRows() ; i++ )
	   {
		// *position(0) = *X_environment(i);
		*position(0) = *X_environment(0) - dx_stroma/2.0 + i*voxel_size + half_voxel_size;

		double E = evaluate_Matrix3( ECM, X_environment , Y_environment, Z_environment , position );
		double BV = normalizer * evaluate_Matrix3( OxygenSourceHD, X_environment , Y_environment, Z_environment ,
	 position ); if( isnan( BV ) ) { BV = 0.0; }

		vector<string> Colors;
		Colors = hematoxylin_and_eosin_stroma_coloring( E , BV );
		Write_SVG_rect( os , *position(0)-half_voxel_size-X_lower , *position(1)-half_voxel_size+top_margin-Y_lower,
		voxel_size , voxel_size , 1 , Colors[0], Colors[0] );

	   }
	  }

	 }
	*/
	os << "  </g>" << std::endl;

	// Now draw vessels

	/*
	 std::vector<std::string> VesselColors = hematoxylin_and_eosin_stroma_coloring( 0,1 );

	 os << " <g id=\"BloodVessels\">" << endl;
	 extern vector<BloodVesselSegment*> BloodVesselSegments;
	 Vector Offset;
	 *Offset(0) = X_lower;
	 *Offset(1) = Y_lower-top_margin;
	*/



	// plot intersecting cells
	os << "  <g id=\"cells\">" << std::endl;
	for (int i = 0; i < total_cell_count; i++)
	{
		cell* pC = e.cast_container<cell_container>().get_at(i); // global_cell_list[i];

		static std::vector<std::string> Colors;
		if (fabs((pC->get_position())[2] - z_slice) < pC->phenotype.geometry.radius())
		{
			double r = pC->phenotype.geometry.radius();
			double rn = pC->phenotype.geometry.nuclear_radius();
			double z = fabs((pC->position())[2] - z_slice);

			Colors = cell_coloring_function(pC);

			os << "   <g id=\"cell" << pC->id << "\" "
			   << "type=\"" << pC->type_name << "\" "; // new April 2022
			if (pC->phenotype.death.dead() == true)
			{
				os << "dead=\"true\" ";
			}
			else
			{
				os << "dead=\"false\" ";
			}
			os << ">" << std::endl;

			// figure out how much of the cell intersects with z = 0

			double plot_radius = sqrt(r * r - z * z);

			Write_SVG_circle(os, (pC->position())[0] - X_lower, (pC->position())[1] - Y_lower, plot_radius, 0.5,
							 Colors[1], Colors[0]);

			// plot the nucleus if it, too intersects z = 0;
			if (fabs(z) < rn && PhysiCell_SVG_options.plot_nuclei == true)
			{
				plot_radius = sqrt(rn * rn - z * z);
				Write_SVG_circle(os, (pC->position())[0] - X_lower, (pC->position())[1] - Y_lower, plot_radius, 0.5,
								 Colors[3], Colors[2]);
			}
			os << "   </g>" << std::endl;
		}
	}
	os << "  </g>" << std::endl;

	// plot intersecting BM points
	/*
	 for( int i=0 ; i < BasementMembraneNodes.size() ; i++ )
	 {
		// vector<string> Colors = false_cell_coloring( pC );
		BasementMembraneNode* pBMN = BasementMembraneNodes[i];
		double thickness =0.1;

		if( fabs( *(pBMN->Position)(2) - z_slice ) < thickness/2.0 )
		{
		 string bm_color ( "rgb(0,0,0)" );
		 double r = thickness/2.0;
		 double z = fabs( *(pBMN->Position)(2) - z_slice) ;

		 os << " <g id=\"BMN" << pBMN->ID << "\">" << std::endl;
		 Write_SVG_circle( os,*(pBMN->Position)(0)-X_lower, *(pBMN->Position)(1)+top_margin-Y_lower, 10*thickness/2.0 ,
	 0.5 , bm_color , bm_color ); os << " </g>" << std::endl;
		}
		// pC = pC->pNextCell;
	 }
	*/

	// end of the <g ID="tissue">
	os << " </g>" << std::endl;

	// draw a scale bar

	double bar_margin = 0.025 * plot_height;
	double bar_height = 0.01 * plot_height;
	double bar_width = PhysiCell_SVG_options.length_bar;
	// double bar_stroke_width = 0.001 * plot_height;

	std::string bar_units = PhysiCell_SVG_options.simulation_space_units;
	// convert from micron to mm
	double temp = bar_width;

	if (temp > 999 && std::strstr(bar_units.c_str(), PhysiCell_SVG_options.mu.c_str()))
	{
		temp /= 1000;
		bar_units = "mm";
	}
	// convert from mm to cm
	if (temp > 9 && std::strcmp(bar_units.c_str(), "mm") == 0)
	{
		temp /= 10;
		bar_units = "cm";
	}

	szString = new char[1024];
	sprintf(szString, "%u %s", (int)round(temp), bar_units.c_str());

	Write_SVG_rect(os, plot_width - bar_margin - bar_width, plot_height + top_margin - bar_margin - bar_height,
				   bar_width, bar_height, 0.002 * plot_height, "rgb(255,255,255)", "rgb(0,0,0)");
	Write_SVG_text(os, szString, plot_width - bar_margin - bar_width + 0.25 * font_size,
				   plot_height + top_margin - bar_margin - bar_height - 0.25 * font_size, font_size,
				   PhysiCell_SVG_options.font_color.c_str(), PhysiCell_SVG_options.font.c_str());

	delete[] szString;

	// plot runtime
	szString = new char[1024];
	RUNTIME_TOC();
	std::string formatted_stopwatch_value = format_stopwatch_value(runtime_stopwatch_value());
	Write_SVG_text(os, formatted_stopwatch_value.c_str(), bar_margin, top_margin + plot_height - bar_margin,
				   0.75 * font_size, PhysiCell_SVG_options.font_color.c_str(), PhysiCell_SVG_options.font.c_str());
	delete[] szString;

	// draw a box around the plot window
	Write_SVG_rect(os, 0, top_margin, plot_width, plot_height, 0.002 * plot_height, "rgb(0,0,0)", "none");

	// close the svg tag, close the file
	Write_SVG_end(os);
	os.close();

	return;
}

void create_plot_legend(std::string filename, cell_coloring_funct_t cell_coloring_function, environment& e)
{
	int number_of_cell_types = e.cell_definitions_count;

	double temp_cell_radius = 25;
	// double temp_cell_volume = 4.1887902047863909846168578443727 * pow(temp_cell_radius, 3.0);

	double relative_padding = 0.15;
	double padding = relative_padding * 2.0 * temp_cell_radius;

	double row_height = 2.0 * temp_cell_radius + 2 * padding;

	double font_size = 0.85 * 2.0 * temp_cell_radius;
	double row_width = 2.0 * temp_cell_radius + 2 * padding + (32 * font_size) + 2 * padding;

	double total_height = number_of_cell_types * row_height;
	double total_width = row_width;


	std::ofstream os(filename, std::ios::out);
	Write_SVG_start(os, total_width, total_height);

	double cursor_x = padding + temp_cell_radius;
	double cursor_y = padding + temp_cell_radius;

	for (int k = 0; k < number_of_cell_types; k++)
	{
		// switch to the cell type
		cell* C = e.cast_container<cell_container>().create_cell(e.cell_definitions[k]);

		// get the colors using the current coloring function
		std::vector<std::string> colors = cell_coloring_function(C);

		// place a big circle with cytoplasm colors
		Write_SVG_circle(os, cursor_x, cursor_y, temp_cell_radius, 1.0, colors[1], colors[0]);
		// place a small circle with nuclear colors
		Write_SVG_circle(os, cursor_x, cursor_y, 0.5 * temp_cell_radius, 1.0, colors[3], colors[2]);

		// place the label

		cursor_x += temp_cell_radius + 2 * padding;
		cursor_y += 0.3 * font_size;

		Write_SVG_text(os, e.cell_definitions[k].name.c_str(), cursor_x, cursor_y, font_size,
					   PhysiCell_SVG_options.font_color.c_str(), PhysiCell_SVG_options.font.c_str());

		// move the cursor down to the next row

		cursor_y -= 0.3 * font_size;
		cursor_y += (2.0 * padding + 2.0 * temp_cell_radius);
		cursor_x = padding + temp_cell_radius;

		e.cast_container<cell_container>().remove_agent(C);
	}

	Write_SVG_end(os);
	os.close();
}

std::vector<std::string> paint_by_number_cell_coloring(cell* pCell)
{
	static std::vector<std::string> colors(0);
	static bool setup_done = false;
	if (setup_done == false)
	{
		colors.push_back("grey"); // default color will be grey

		colors.push_back("red");
		colors.push_back("yellow");
		colors.push_back("green");
		colors.push_back("blue");

		colors.push_back("magenta");
		colors.push_back("orange");
		colors.push_back("lime");
		colors.push_back("cyan");

		colors.push_back("hotpink");
		colors.push_back("peachpuff");
		colors.push_back("darkseagreen");
		colors.push_back("lightskyblue");

		setup_done = true;
	}

	// start all black

	std::vector<std::string> output = { "black", "black", "black", "black" };

	// paint by number -- by cell type

	std::string interior_color = "white";
	if (pCell->type < 13)
	{
		interior_color = colors[pCell->type];
	}

	output[0] = interior_color; // set cytoplasm color

	/*
	if( pCell->phenotype.death.dead == false ) // if live, color nucleus same color
	{
		output[2] = interior_color;
		output[3] = interior_color;
	}
	else
	{
		// apoptotic cells will retain a black nucleus
		// if necrotic, color the nucleus brown
		if( pCell->phenotype.cycle.current_phase().code == constants::necrotic_swelling ||
			pCell->phenotype.cycle.current_phase().code == constants::necrotic_lysed ||
			pCell->phenotype.cycle.current_phase().code == constants::necrotic )
		{
			output[2] = "rgb(139,69,19)";
			output[3] = "rgb(139,69,19)";
		}



	}
	*/

	// new March 2023 (for better compatibility with studio)

	// if dead, use live color for the outline
	// if( pCell->phenotype.death.dead == true )
	// { output[1] = interior_color; }

	// necrotic cells are brown
	if (pCell->phenotype.cycle.current_phase().code == constants::necrotic_swelling
		|| pCell->phenotype.cycle.current_phase().code == constants::necrotic_lysed
		|| pCell->phenotype.cycle.current_phase().code == constants::necrotic)
	{
		interior_color = "saddlebrown";
	}
	// apoptotic cells are white
	if (pCell->phenotype.cycle.current_phase().code == constants::apoptotic)
	{
		interior_color = "black";
	}

	output[0] = interior_color; // set cytoplasm color
	output[2] = interior_color; // set cytoplasm color
	output[3] = interior_color; // set cytoplasm color

	output[1] = "black";

	return output;
}

} // namespace physicell
