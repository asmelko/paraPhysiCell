
#include "custom_cell_data.h"

#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

using namespace biofvm;
using namespace physicell;

Variable::Variable()
{
	name = "unnamed";
	units = "dimensionless";
	value = 0.0;
	conserved_quantity = false;
	return;
}

std::ostream& operator<<(std::ostream& os, const Variable& v)
{
	os << v.name << ": " << v.value << " " << v.units;
	return os;
}


Vector_Variable::Vector_Variable()
{
	name = "unnamed";
	units = "dimensionless";
	value.resize(3, 0.0);
	conserved_quantity = false;
	return;
}

std::ostream& operator<<(std::ostream& os, const Vector_Variable& v)
{
	os << v.name << ": ";
	if (v.value.size() == 0)
	{
		os << "[empty]";
		return os;
	}
	/*
		if( v.value.size() == 1 )
		{ os << v.value[0] << " (" << v.units << ")"; return os;  }
	*/
	for (std::size_t i = 0; i < v.value.size() - 1; i++)
	{
		os << v.value[i] << ",";
	}
	os << v.value[v.value.size() - 1] << " (" << v.units << ")";
	return os;
}

Custom_Cell_Data::Custom_Cell_Data(const Custom_Cell_Data& ccd)
{
	//	std::cout << __FUNCTION__ << "(copy)" << std::endl;
	variables = ccd.variables;
	vector_variables = ccd.vector_variables;

	name_to_index_map = ccd.name_to_index_map;

	return;
}

index_t Custom_Cell_Data::add_variable(Variable& v)
{
	index_t n = variables.size();
	variables.push_back(v);
	name_to_index_map[v.name] = n;
	return n;
}

index_t Custom_Cell_Data::add_variable(std::string name, std::string units, real_t value)
{
	index_t n = variables.size();
	variables.resize(n + 1);
	variables[n].name = name;
	variables[n].units = units;
	variables[n].value = value;
	name_to_index_map[name] = n;
	return n;
}

index_t Custom_Cell_Data::add_variable(std::string name, real_t value)
{
	index_t n = variables.size();
	variables.resize(n + 1);
	variables[n].name = name;
	variables[n].units = "dimensionless";
	variables[n].value = value;
	name_to_index_map[name] = n;
	return n;
}

index_t Custom_Cell_Data::add_vector_variable(Vector_Variable& v)
{
	index_t n = vector_variables.size();
	vector_variables.push_back(v);
	//	vector_name_to_index_map[ v.name ] = n;
	return n;
}

index_t Custom_Cell_Data::add_vector_variable(std::string name, std::string units, std::vector<real_t>& value)
{
	index_t n = vector_variables.size();
	vector_variables.resize(n + 1);
	vector_variables[n].name = name;
	vector_variables[n].units = units;
	vector_variables[n].value = value;
	//	vector_name_to_index_map[ name ] = n;
	return n;
}

index_t Custom_Cell_Data::add_vector_variable(std::string name, std::vector<real_t>& value)
{
	index_t n = vector_variables.size();
	vector_variables.resize(n + 1);
	vector_variables[n].name = name;
	vector_variables[n].units = "dimensionless";
	vector_variables[n].value = value;
	//	vector_name_to_index_map[ name ] = n;
	return n;
}

index_t Custom_Cell_Data::find_variable_index(std::string name)
{
	// this should return -1 if not found, not zero
	auto out = name_to_index_map.find(name);
	if (out != name_to_index_map.end())
	{
		return out->second;
	}
	return -1;
}

/*
index_t Custom_Cell_Data::find_vector_variable_index( std::string name )
{
	return vector_name_to_index_map[ name ];
}
*/

index_t Custom_Cell_Data::find_vector_variable_index(std::string name)
{
	std::size_t n = 0;
	while (n < vector_variables.size())
	{
		if (std::strcmp(vector_variables[n].name.c_str(), name.c_str()) == 0)
		{
			return n;
		}
		n++;
	}

	return -1;
}

real_t& Custom_Cell_Data::operator[](index_t i) { return variables[i].value; }

real_t& Custom_Cell_Data::operator[](std::string name) { return variables[name_to_index_map[name]].value; }

std::ostream& operator<<(std::ostream& os, const Custom_Cell_Data& ccd)
{
	os << "Custom data (scalar): " << std::endl;
	for (std::size_t i = 0; i < ccd.variables.size(); i++)
	{
		os << i << ": " << ccd.variables[i] << std::endl;
	}

	os << "Custom data (vector): " << std::endl;
	for (std::size_t i = 0; i < ccd.vector_variables.size(); i++)
	{
		os << i << ": " << ccd.vector_variables[i] << std::endl;
	}

	return os;
}
