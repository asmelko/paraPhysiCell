#pragma once

#include <pugixml.hpp>

#include <BioFVM/types.h>

namespace physicell {

// find the first <find_me> child in <parent_node>
pugi::xml_node xml_find_node(pugi::xml_node& parent_node, std::string find_me); // done

// get the std:string in <parent_node> <find_me>string_value</find_me> </parent_node>
std::string xml_get_string_value(pugi::xml_node& parent_node, std::string find_me); // done

// get the biofvm::real_t value stored in <parent_node> <find_me>double_value</find_me> </parent_node>
biofvm::real_t xml_get_double_value(pugi::xml_node& parent_node, std::string find_me); // done

// get the integer value in <parent_node> <find_me>int_value</find_me> </parent_node>
biofvm::index_t xml_get_int_value(pugi::xml_node& parent_node, std::string find_me); // done

// get the Boolean value in <parent_node> <find_me>int_value</find_me> </parent_node>
bool xml_get_bool_value(pugi::xml_node& parent_node, std::string find_me); // done


// get the name of the element in <my_node> (the name would be my_node)
std::string xml_get_my_name(pugi::xml_node node);


bool xml_get_my_bool_value(pugi::xml_node node);
biofvm::index_t xml_get_my_int_value(pugi::xml_node node);
biofvm::real_t xml_get_my_double_value(pugi::xml_node node);
std::string xml_get_my_string_value(pugi::xml_node node);



// get the string attribute named "attribute" in the first std:string in <parent_node> <find_me>string_value</find_me>
// </parent_node>
std::string get_string_attribute_value(pugi::xml_node& parent_node, std::string find_me, std::string attribute);

// get the biofvm::index_t attribute named "attribute" in the first std:string in <parent_node>
// <find_me>string_value</find_me>
// </parent_node>
biofvm::index_t get_int_attribute_value(pugi::xml_node& parent_node, std::string find_me, std::string attribute);

// get the biofvm::real_t attribute named "attribute" in the first std:string in <parent_node>
// <find_me>string_value</find_me>
// </parent_node>
biofvm::real_t get_double_attribute_value(pugi::xml_node& parent_node, std::string find_me, std::string attribute);

}; // namespace physicell
