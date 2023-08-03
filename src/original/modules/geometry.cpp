#include "geometry.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "../signal_behavior.h"
#include "pugixml_helper.h"
#include "types.h"
#include "vector_utils.h"

using namespace biofvm;

namespace physicell {

std::vector<std::string> split_csv_labels(std::string labels_line)
{
	std::vector<std::string> label_tokens;
	std::string s;

	std::stringstream stream(labels_line);
	while (std::getline(stream, s, ','))
	{
		label_tokens.push_back(s);
	}

	return label_tokens;
}

void load_cells_csv_v1(std::string filename, environment& e)
{
	std::ifstream file(filename, std::ios::in);
	if (!file)
	{
		std::cout << "Error: " << filename << " not found during cell loading. Quitting." << std::endl;
		exit(-1);
	}
	else
	{
		std::cout << "Loading cells from simple (v1) CSV file " << filename << " ... " << std::endl;
	}

	std::string line;
	while (std::getline(file, line))
	{
		std::vector<real_t> data = csv_to_vector(line.c_str());

		if (data.size() != 4)
		{
			std::cout << "Error! Importing cells from a CSV file expects each row to be x,y,z,typeID." << std::endl;
			exit(-1);
		}

		point_t<real_t, 3> position = { data[0], data[1], data[2] };

		int my_type = (int)data[3];
		cell_definition* pCD = e.find_cell_definition(my_type);
		if (pCD != NULL)
		{
			std::cout << "Creating " << pCD->name << " (type=" << pCD->type << ") at " << position << std::endl;
			cell* pCell = e.cast_container<cell_container>().create_cell(*pCD);
			pCell->assign_position(position);
		}
		else
		{
			std::cout << "Warning! No cell definition found for index " << my_type << "!" << std::endl
					  << "\tIgnoring cell in " << filename << " at position " << position << std::endl;
		}
	}

	file.close();
}


cell* process_csv_v2_line(std::string line, std::vector<std::string> labels, environment& e)
{
	// split the line into tokens
	std::vector<std::string> tokens;

	std::stringstream stream(line);
	std::string s;
	while (std::getline(stream, s, ','))
	{
		tokens.push_back(s);
	}

	// get the cell position
	point_t<real_t, 3> position;
	char* pTemp;
	for (int i = 0; i < 3; i++)
	{
		position[i] = strtod(tokens[i].c_str(), &pTemp);
	}

	// the cell type
	std::string celltype = tokens[3];
	cell_definition* pCD = e.find_cell_definition(celltype);
	if (pCD == NULL)
	{
		std::cout << "Warning! CSV file requests creating cell type " << celltype << std::endl
				  << "\tat " << position << "but I don't recognize that type. Skipping cell!" << std::endl
				  << std::endl;
		return NULL;
	}

	// create the cell IF the definition was found
	std::cout << "Creating " << pCD->name << " (type=" << pCD->type << ") at " << position << std::endl;

	cell* pCell = e.cast_container<cell_container>().create_cell(*pCD);
	pCell->assign_position(position);

	// now write any extra data

	for (std::size_t k = 4; k < tokens.size(); k++)
	{
		double dval = strtod(tokens[k].c_str(), &pTemp);
		bool processed = false;
		bool skip = false;


		// if the string is empty, skip
		if (tokens[k].size() == 0)
		{
			skip = true;
		}
		else
		{
			char c = tokens[k].c_str()[0];
			// use 's' or 'S' to skip the entry
			if (c == 's' || c == 'S')
			{
				skip = true;
			}
		}

		// special cases:

		// volume
		if (labels[k] == "volume" && skip == false)
		{
			pCell->set_total_volume(dval);
			processed = true;
		}

		// check behavior dictionary

		if (processed == false && skip == false)
		{
			// if the behavior is found in the dictionary, process it
			if (find_behavior_index(labels[k]) > -1)
			{
				set_single_behavior(pCell, labels[k], dval, e);
				processed = true;
			}
		}

		// warning message for any unprocessed variables
		if (processed == false && skip == false)
		{
			std::cout << "\tWarning: I don't know how to process " << labels[k] << " so I skipped it." << std::endl;
		}
		// give a notation for any intentinoally skipped variables
		if (skip == true)
		{
			std::cout << "\tNote: Skipping " << labels[k] << " for this cell." << std::endl;
		}
	}

	return pCell;
}

void load_cells_csv_v2(std::string filename, environment& e)
{
	// open file
	std::ifstream file(filename, std::ios::in);
	if (!file)
	{
		std::cout << "Error: " << filename << " not found during cell loading. Quitting." << std::endl;
		exit(-1);
	}
	else
	{
		std::cout << "Loading cells from detailed (v2) CSV file " << filename << " ... " << std::endl;
	}

	// get the first line (labels)

	std::string line;
	std::getline(file, line);

	// tokenize the labels

	std::vector<std::string> labels = split_csv_labels(line);

	// process all remaining lines

	while (std::getline(file, line))
	{
		process_csv_v2_line(line, labels, e);
	}

	// close the file

	file.close();
	std::cout << "Done! " << std::endl << std::endl;

	return;
}

void load_cells_csv(std::string filename, environment& e)
{
	// open file
	std::ifstream file(filename, std::ios::in);
	if (!file)
	{
		std::cout << "Error: " << filename << " not found during cell loading. Quitting." << std::endl;
		exit(-1);
	}

	// determine version
	std::string line;
	std::getline(file, line);
	char c = line.c_str()[0];

	file.close();

	if (c == 'X' || c == 'x')
	{
		// v2
		return load_cells_csv_v2(filename, e);
	}
	else
	{
		// v1
		return load_cells_csv_v1(filename, e);
	}

	return;
}


bool load_cells_from_pugixml(const pugi::xml_node& root, environment& e)
{
	pugi::xml_node node = root.child("initial_conditions");
	if (!node)
	{
		std::cout << "Warning: XML-based cell positions has wrong formating. Ignoring!" << std::endl;
		return false;
	}

	node = node.child("cell_positions");
	if (!node)
	{
		std::cout << "Warning: XML-based cell positions has wrong formating. Ignoring!" << std::endl;
		return false;
	}

	// enabled?
	if (node.attribute("enabled").as_bool() == false)
	{
		return false;
	}

	// get filename

	std::string folder = xml_get_string_value(node, "folder");
	std::string filename = xml_get_string_value(node, "filename");
	std::string input_filename = folder + "/" + filename;

	std::string filetype = node.attribute("type").value();

	// what kind?
	if (filetype == "csv" || filetype == "CSV")
	{
		std::cout << "Loading cells from CSV file " << input_filename << " ... " << std::endl;
		load_cells_csv(input_filename, e);
		system("sleep 1");
		return true;
	}
	if (filetype == "matlab" || filetype == "mat" || filetype == "MAT")
	{
		std::cout << "Error: Load cell positions from matlab not yet supported. Try CSV." << std::endl;
		exit(-1);
		std::cout << "Loading cells from matlab file " << input_filename << " ... " << std::endl;
		return false;
	}
	if (filetype == "scene")
	{
		std::cout << "Error: load cell positions from scene not yet supported. Try CSV." << std::endl;
		exit(-1);
		std::cout << "Loading cells from scene file " << input_filename << " ... " << std::endl;
		return false;
	}
	if (filetype == "physicell" || filetype == "PhysiCell")
	{
		std::cout << "Error: load cell positions from PhysiCell snapshot not yet supported. Try CSV." << std::endl;
		exit(-1);
		std::cout << "Loading cells from PhysiCell file " << input_filename << " ... " << std::endl;
		return false;
	}

	return false;
}

} // namespace physicell
