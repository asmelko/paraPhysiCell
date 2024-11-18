#include "settings.h"

#include <iostream>

#include "../../random.h"
#include "pugixml_helper.h"

using namespace biofvm;

namespace physicell {

PhysiCell_Settings::PhysiCell_Settings()
{
	// units
	time_units = "min";
	space_units = "micron";

	// save options
	folder = ".";
	max_time = 60 * 24 * 45;

	full_save_interval = 60;
	enable_full_saves = true;
	enable_legacy_saves = false;

	SVG_save_interval = 60;
	enable_SVG_saves = true;
	enable_substrate_plot = false;
	substrate_to_monitor = "oxygen";
	limits_substrate_plot = false;
	min_concentration = -1.0;
	max_concentration = -1.0;
	svg_substrate_colormap = "YlOrRd";

	intracellular_save_interval = 60;
	enable_intracellular_saves = false;

	// parallel options

	omp_num_threads = 4;

	rules_enabled = false;
	rules_protocol = "Cell Behavior Hypothesis Grammar (CBHG)";
	rules_protocol_version = "1.0";

	return;
}

void PhysiCell_Settings::read_from_pugixml(pugi::xml_node& physicell_config_root)
{
	pugi::xml_node node;

	// overall options

	node = xml_find_node(physicell_config_root, "overall");

	max_time = xml_get_double_value(node, "max_time");
	time_units = xml_get_string_value(node, "time_units");
	space_units = xml_get_string_value(node, "space_units");

	node = node.parent();

	// save options

	node = xml_find_node(physicell_config_root, "save");

	folder = xml_get_string_value(node, "folder");

	node = xml_find_node(node, "full_data");
	full_save_interval = xml_get_double_value(node, "interval");
	enable_full_saves = xml_get_bool_value(node, "enable");
	node = node.parent();

	node = xml_find_node(node, "SVG");
	SVG_save_interval = xml_get_double_value(node, "interval");
	enable_SVG_saves = xml_get_bool_value(node, "enable");

	pugi::xml_node node_plot_substrate;
	node_plot_substrate = xml_find_node(node, "plot_substrate");
	enable_substrate_plot = node_plot_substrate.attribute("enabled").as_bool();
	limits_substrate_plot = node_plot_substrate.attribute("limits").as_bool();

	if (enable_substrate_plot)
	{
		substrate_to_monitor = xml_get_string_value(node_plot_substrate, "substrate");
		if (limits_substrate_plot)
		{
			min_concentration = xml_get_double_value(node_plot_substrate, "min_conc");
			max_concentration = xml_get_double_value(node_plot_substrate, "max_conc");
		}
		pugi::xml_node colormap_node = xml_find_node(node_plot_substrate, "colormap");
		if (colormap_node)
		{
			svg_substrate_colormap = xml_get_my_string_value(colormap_node);
		}
	}

	node = node.parent();

	node = xml_find_node(node, "intracellular_data");
	intracellular_save_interval = xml_get_double_value(node, "interval");
	enable_intracellular_saves = xml_get_bool_value(node, "enable");
	node = node.parent();

	node = xml_find_node(node, "legacy_data");
	enable_legacy_saves = xml_get_bool_value(node, "enable");
	node = node.parent();

	// parallel options

	node = xml_find_node(physicell_config_root, "parallel");
	omp_num_threads = xml_get_int_value(node, "omp_num_threads");

	node = node.parent();

	// legacy and other options

	pugi::xml_node node_options;

	node_options = xml_find_node(physicell_config_root, "options");
	if (node_options)
	{
		bool settings;

		// look for legacy_random_points_on_sphere_in_divide
		settings = xml_get_bool_value(node_options, "legacy_random_points_on_sphere_in_divide");
		if (settings)
		{
			std::cout << "setting legacy unif" << std::endl;
			throw std::runtime_error("Error: legacy_random_points_on_sphere_in_divide is not supported.");
			// extern std::vector<double> (*cell_division_orientation)(void);
			// cell_division_orientation = LegacyRandomOnUnitSphere;
		}

		settings = xml_get_bool_value(node_options, "disable_automated_spring_adhesions");
		if (settings)
		{
			std::cout << "Disabling automated spring adhesions and detachments!" << std::endl;
			disable_automated_spring_adhesions = true;
		}

		pugi::xml_node random_seed_node = xml_find_node(node_options, "random_seed");
		std::string random_seed = ""; // default is system clock, even if this element is not present
		if (random_seed_node)
		{
			random_seed = xml_get_my_string_value(random_seed_node);
		}

		if (random_seed == "" || random_seed == "random" || random_seed == "system_clock")
		{
			std::cout << "Using system clock for random seed" << std::endl;
			random().set_random_seed();
		}
		else
		{
			int seed;
			try
			{
				seed = std::stoi(random_seed);
			}
			catch (const std::exception& e)
			{
				std::cout << "ERROR: " << random_seed
						  << " is not a valid random seed. It must be an integer. Fix this within <options>."
						  << std::endl;
				exit(-1);
			}
			random().set_seed(seed);
		}

		// other options can go here, eventually
	}
}

template <class T>
Parameter<T>::Parameter()
{
	name = "unnamed";
	units = "none";
	/*
		T* pT;
		pT = new T;
		value = *pT;
	*/
	value = (T)0;
	//	value = 1-1;
	return;
}

template <>
Parameter<std::string>::Parameter()
{
	name = "unnamed";
	units = "none";
	value = "none";
	return;
}

template <class T>
Parameter<T>::Parameter(std::string my_name)
{
	name = my_name;
	units = "dimensionless";
	/*
		T* pT;
		pT = new T;
		value = *pT;
	*/
	value = (T)0;
	return;
}

template <>
Parameter<std::string>::Parameter(std::string my_name)
{
	name = my_name;
	units = "none";
	value = "none";
	return;
}

template <class T>
void Parameter<T>::operator=(T& rhs)
{
	value = rhs;
	return;
}

template <class T>
void Parameter<T>::operator=(T rhs)
{
	value = rhs;
	return;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Parameter<T>& param)
{
	os << param.name << ": " << param.value << " [" << param.units << "]";
	return os;
}

template <class T>
index_t Parameters<T>::size() const
{
	return parameters.size();
}

template <class T>
T& Parameters<T>::operator()(index_t i)
{
	return parameters[i].value;
}

template <class T>
T& Parameters<T>::operator()(std::string str)
{
	if (name_to_index_map.find(str) == name_to_index_map.end())
	{
		std::cerr << "ERROR : Unknown parameter " << str << " ! Quitting." << std::endl;
		exit(-1);
	}
	return parameters[name_to_index_map[str]].value;
}

template <class T>
Parameter<T>& Parameters<T>::operator[](index_t i)
{
	return parameters[i];
}

template <class T>
Parameter<T>& Parameters<T>::operator[](std::string str)
{
	if (name_to_index_map.find(str) == name_to_index_map.end())
	{
		std::cerr << "ERROR : Unknown parameter " << str << " ! Quitting." << std::endl;
		exit(-1);
	}
	return parameters[name_to_index_map[str]];
}


template <class T>
index_t Parameters<T>::find_index(std::string search_name)
{
	auto out = name_to_index_map.find(search_name);
	if (out != name_to_index_map.end())
	{
		return out->second;
	}
	return -1;
	// return name_to_index_map[ search_name ];
}


template <class T>
std::ostream& operator<<(std::ostream& os, const Parameters<T>& params)
{
	for (std::size_t i = 0; i < params.parameters.size(); i++)
	{
		os << params.parameters[i] << std::endl;
	}
	return os;
}

template <class T>
Parameters<T>::Parameters()
{
	parameters.resize(0);
	name_to_index_map.clear();
	return;
}

template <class T>
void Parameters<T>::add_parameter(std::string my_name)
{
	// this function is not currently (2024-06-03) called in the code, so these defaults largely do not matter; very
	// unlikely others are directly calling this function, let alone this implementation
	T my_value = T(); // for {int, double, bool, string} this will be {0, 0.0, false, ""} (this would technically change
					  // the behavior for strings since it is hardcoded above to default to "none", but nobody should
					  // rely on the default value of a string being "none")
	return add_parameter(my_name, my_value);
}

template <class T>
void Parameters<T>::add_parameter(std::string my_name, T my_value)
{
	// this function is not currently (2024-06-03) called in the code, so these defaults largely do not matter; very
	// unlikely others are directly calling this function, let alone this implementation
	std::string my_units =
		"dimensionless"; // technically this would change the behavior for strings since it is hardcoded above to
						 // default to "none", but nobody should be using units on strings; also, if the xml does not
						 // have units, then "dimensionless" is used even for strings
	return add_parameter(my_name, my_value, my_units);
}

template <class T>
void Parameters<T>::add_parameter(std::string my_name, T my_value, std::string my_units)
{
	assert_not_exists(my_name);

	Parameter<T> new_param;
	new_param.name = my_name;
	new_param.value = my_value;
	new_param.units = my_units;

	index_t n = parameters.size();

	parameters.push_back(std::move(new_param));

	name_to_index_map[my_name] = n;
	return;
}

template <class T>
void Parameters<T>::add_parameter(Parameter<T> param)
{
	assert_not_exists(param.name);

	int n = parameters.size();
	parameters.push_back(param);
	name_to_index_map[param.name] = n;
	return;
}

template <class T>
void Parameters<T>::assert_not_exists(std::string search_name)
{
	if (find_index(search_name) == -1)
	{
		return;
	}

	std::cout << "ERROR: Parameter " << search_name
			  << " already exists. Make sure all parameters (of a given type) have unique names." << std::endl;
	exit(-1);
}

std::ostream& operator<<(std::ostream& os, const User_Parameters up)
{
	os << "Bool parameters:: " << std::endl << up.bools << std::endl;
	os << "Int parameters:: " << std::endl << up.ints << std::endl;
	os << "Double parameters:: " << std::endl << up.doubles << std::endl;
	os << "String parameters:: " << std::endl << up.strings << std::endl;
	return os;
}

void User_Parameters::read_from_pugixml(pugi::xml_node& parent_node)
{
	pugi::xml_node node = xml_find_node(parent_node, "user_parameters");

	pugi::xml_node node1 = node.first_child();
	while (node1)
	{
		std::string name = xml_get_my_name(node1);
		std::string units = node1.attribute("units").value();
		if (units == "")
		{
			units = "dimensionless";
		}

		std::string type = node1.attribute("type").value();

		if (type == "bool")
		{
			bool value = xml_get_my_bool_value(node1);
			bools.add_parameter(name, value, units);
		}
		else if (type == "int")
		{
			int value = xml_get_my_int_value(node1);
			ints.add_parameter(name, value, units);
		}
		else if (type == "double")
		{
			double value = xml_get_my_double_value(node1);
			doubles.add_parameter(name, value, units);
		}
		else if (type == "string")
		{
			std::string value = xml_get_my_string_value(node1);
			strings.add_parameter(name, value, units);
		}
		else // default if no type specified
		{
			double value = xml_get_my_double_value(node1);
			doubles.add_parameter(name, value, units);
		}

		node1 = node1.next_sibling();
	}

	std::cout << "User parameters in XML config file: " << std::endl;
	std::cout << *this << std::endl;

	return;
}

// need this so that the template gets filled and compiled prior to linking
template class Parameter<bool>;
template class Parameter<index_t>;
template class Parameter<real_t>;
template class Parameter<std::string>;

template class Parameters<bool>;
template class Parameters<index_t>;
template class Parameters<real_t>;
template class Parameters<std::string>;

} // namespace physicell
