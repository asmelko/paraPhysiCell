#pragma once

#include <pugixml.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include <BioFVM/types.h>

namespace physicell {

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
	void operator=(Parameter& p);
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

	biofvm::index_t size(void) const;
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

	void read_from_pugixml(pugi::xml_node parent_node);
};

} // namespace physicell
