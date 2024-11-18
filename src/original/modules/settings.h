#pragma once

#include <pugixml.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include <BioFVM/types.h>

namespace physicell {

class PhysiCell_Settings
{
private:
public:
	// overall
	biofvm::real_t max_time = 60 * 24 * 45;

	// units
	std::string time_units = "min";
	std::string space_units = "micron";

	// parallel options
	biofvm::index_t omp_num_threads = 2;

	// save options
	std::string folder = ".";

	biofvm::real_t full_save_interval = 60;
	bool enable_full_saves = true;
	bool enable_legacy_saves = false;

	bool disable_automated_spring_adhesions = false;

	biofvm::real_t SVG_save_interval = 60;
	bool enable_SVG_saves = true;

	bool enable_substrate_plot = false;
	std::string substrate_to_monitor = "oxygen";
	bool limits_substrate_plot = false;
	biofvm::real_t min_concentration = -1.0;
	biofvm::real_t max_concentration = -1.0;
	std::string svg_substrate_colormap = "YlOrRd";

	biofvm::real_t intracellular_save_interval = 60;
	bool enable_intracellular_saves = false;

	// cell rules option
	bool rules_enabled = false;
	std::string rules_protocol = "Cell Behavior Hypothesis Grammar (CBHG)";
	std::string rules_protocol_version = "1.0";

	PhysiCell_Settings();

	void read_from_pugixml(pugi::xml_node& physicell_config_root);
};

template <class T>
class Parameter
{
private:
	template <class Y>
	friend std::ostream& operator<<(std::ostream& os, const Parameter<Y>& param);

public:
	std::string name;
	std::string units;
	T value;

	Parameter();
	Parameter(std::string my_name);

	void operator=(T& rhs);
	void operator=(T rhs);
};

template <class T>
class Parameters
{
private:
	std::unordered_map<std::string, biofvm::index_t> name_to_index_map;

	template <class Y>
	friend std::ostream& operator<<(std::ostream& os, const Parameters<Y>& params);

public:
	Parameters();

	std::vector<Parameter<T>> parameters;

	void add_parameter(std::string my_name);
	void add_parameter(std::string my_name, T my_value);
	//	void add_parameter( std::string my_name , T my_value );
	void add_parameter(std::string my_name, T my_value, std::string my_units);
	//	void add_parameter( std::string my_name , T my_value , std::string my_units );

	void add_parameter(Parameter<T> param);

	biofvm::index_t find_index(std::string search_name);

	// these access the values
	T& operator()(biofvm::index_t i);
	T& operator()(std::string str);

	// these access the full, raw parameters
	Parameter<T>& operator[](biofvm::index_t i);
	Parameter<T>& operator[](std::string str);

	biofvm::index_t size() const;

	void assert_not_exists(std::string search_name);
};

class User_Parameters
{
private:
	friend std::ostream& operator<<(std::ostream& os, const User_Parameters up);

public:
	Parameters<bool> bools;
	Parameters<biofvm::index_t> ints;
	Parameters<biofvm::real_t> doubles;
	Parameters<std::string> strings;

	void read_from_pugixml(pugi::xml_node& parent_node);
};

} // namespace physicell
