#include "user_parameters.h"

#include <iostream>

#include "pugixml_helper.h"

using namespace biofvm;

namespace physicell {

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
void Parameter<T>::operator=(Parameter& p)
{
	name = p.name;
	units = p.units;
	value = p.value;
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
	Parameter<T>* pNew;
	pNew = new Parameter<T>;
	pNew->name = my_name;

	index_t n = parameters.size();

	parameters.push_back(*pNew);

	name_to_index_map[my_name] = n;
	return;
}

template <class T>
void Parameters<T>::add_parameter(std::string my_name, T my_value)
{
	Parameter<T>* pNew;
	pNew = new Parameter<T>;
	pNew->name = my_name;
	pNew->value = my_value;

	index_t n = parameters.size();

	parameters.push_back(*pNew);

	name_to_index_map[my_name] = n;
	return;
}
/*
template <class T>
void Parameters<T>::add_parameter( std::string my_name , T my_value )
{
	Parameter<T>* pNew;
	pNew = new Parameter<T> ;
	pNew->name = my_name ;
	pNew->value = my_value;

	index_t n = parameters.size();

	parameters.push_back( *pNew );

	name_to_index_map[ my_name ] = n;
	return;
}
*/

template <class T>
void Parameters<T>::add_parameter(std::string my_name, T my_value, std::string my_units)
{
	Parameter<T>* pNew;
	pNew = new Parameter<T>;
	pNew->name = my_name;
	pNew->value = my_value;
	pNew->units = my_units;

	index_t n = parameters.size();

	parameters.push_back(*pNew);

	name_to_index_map[my_name] = n;
	return;
}

/*
template <class T>
void Parameters<T>::add_parameter( std::string my_name , T my_value , std::string my_units )
{
	Parameter<T>* pNew;
	pNew = new Parameter<T> ;
	pNew->name = my_name ;
	pNew->value = my_value;
	pNew->units = my_units;

	index_t n = parameters.size();

	parameters.push_back( *pNew );

	name_to_index_map[ my_name ] = n;
	return;
}
*/

template <class T>
void Parameters<T>::add_parameter(Parameter<T> param)
{
	index_t n = parameters.size();
	parameters.push_back(param);
	name_to_index_map[param.name] = n;
	return;
}

std::ostream& operator<<(std::ostream& os, const User_Parameters up)
{
	os << "Bool parameters:: " << std::endl << up.bools << std::endl;
	os << "Int parameters:: " << std::endl << up.ints << std::endl;
	os << "Double parameters:: " << std::endl << up.doubles << std::endl;
	os << "String parameters:: " << std::endl << up.strings << std::endl;
	return os;
}

void User_Parameters::read_from_pugixml(pugi::xml_node parent_node)
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

		bool done = false;
		if (type == "bool" && done == false)
		{
			bool value = xml_get_my_bool_value(node1);
			bools.add_parameter(name, value, units);
			done = true;
		}

		if (type == "int" && done == false)
		{
			index_t value = xml_get_my_int_value(node1);
			ints.add_parameter(name, value, units);
			done = true;
		}

		if (type == "double" && done == false)
		{
			real_t value = xml_get_my_double_value(node1);
			doubles.add_parameter(name, value, units);
			done = true;
		}

		if (done == false && type == "string")
		{
			std::string value = xml_get_my_string_value(node1);
			strings.add_parameter(name, value, units);
			done = true;
		}

		/* default if no type specified: */
		if (done == false)
		{
			real_t value = xml_get_my_double_value(node1);
			doubles.add_parameter(name, value, units);
			done = true;
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