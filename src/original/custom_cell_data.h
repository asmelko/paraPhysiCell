#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <BioFVM/types.h>

namespace physicell {

class Variable
{
public:
	std::string name;
	biofvm::real_t value;
	std::string units;
	bool conserved_quantity;

	Variable();
};

class Vector_Variable
{
public:
	std::string name;
	std::vector<biofvm::real_t> value;
	std::string units;
	bool conserved_quantity;

	Vector_Variable();
};

class Custom_Cell_Data
{
private:
	std::unordered_map<std::string, biofvm::index_t> name_to_index_map;
	//	std::unordered_map<std::string,biofvm::index_t> vector_name_to_index_map;

	friend std::ostream& operator<<(std::ostream& os, const Custom_Cell_Data& ccd); // done
public:
	std::vector<Variable> variables;
	std::vector<Vector_Variable> vector_variables;

	biofvm::index_t add_variable(Variable& v);												 // done
	biofvm::index_t add_variable(std::string name, std::string units, biofvm::real_t value); // done
	biofvm::index_t add_variable(std::string name, biofvm::real_t value);					 // done

	biofvm::index_t add_vector_variable(Vector_Variable& v); // done
	biofvm::index_t add_vector_variable(std::string name, std::string units,
										std::vector<biofvm::real_t>& value);				   // done
	biofvm::index_t add_vector_variable(std::string name, std::vector<biofvm::real_t>& value); // done

	biofvm::index_t find_variable_index(std::string name);		  // done
	biofvm::index_t find_vector_variable_index(std::string name); // done

	// these access the scalar variables
	biofvm::real_t& operator[](biofvm::index_t i); // done
	biofvm::real_t& operator[](std::string name);  // done


	Custom_Cell_Data() = default; // done
	Custom_Cell_Data(const Custom_Cell_Data& ccd);
};

}; // namespace physicell
