#include "pugixml_helper.h"

namespace physicell {


// find the first <find_me> child in <parent_node>
pugi::xml_node xml_find_node(pugi::xml_node& parent_node, std::string find_me)
{
	return parent_node.child(find_me.c_str());
}

// get the std:string in <parent_node> <find_me>string_value</find_me> </parent_node>
std::string xml_get_string_value(pugi::xml_node& parent_node, std::string find_me)
{
	return parent_node.child(find_me.c_str()).text().get();
}


// get the double value stored in <parent_node> <find_me>double_value</find_me> </parent_node>
double xml_get_double_value(pugi::xml_node& parent_node, std::string find_me)
{
	// return strtod( parent_node.child( find_me.c_str() ).text().get() , NULL ); // classic

	return parent_node.child(find_me.c_str()).text().as_double(); // using pugixml conversion
}

// get the integer value in <parent_node> <find_me>int_value</find_me> </parent_node>
int xml_get_int_value(pugi::xml_node& parent_node, std::string find_me)
{
	//	return atoi( parent_node.child( find_me.c_str() ).text().get() ); // classic

	return parent_node.child(find_me.c_str()).text().as_int(); // using pugixml conversion
}

// get the Boolean value in <parent_node> <find_me>int_value</find_me> </parent_node>
bool xml_get_bool_value(pugi::xml_node& parent_node, std::string find_me)
{
	//	return (bool) atoi( parent_node.child( find_me.c_str() ).text().get() ); // classic (untested)

	return parent_node.child(find_me.c_str()).text().as_bool(); // using pugixml conversion
}

// get the name of the element in <my_node> (the name would be my_node)
std::string xml_get_my_name(pugi::xml_node node) { return node.name(); }

bool xml_get_my_bool_value(pugi::xml_node node) { return node.text().as_bool(); }

int xml_get_my_int_value(pugi::xml_node node) { return node.text().as_int(); }

double xml_get_my_double_value(pugi::xml_node node) { return node.text().as_double(); }

std::string xml_get_my_string_value(pugi::xml_node node) { return node.text().get(); }

}; // namespace physicell
