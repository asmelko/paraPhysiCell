#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include <BioFVM/types.h>

#include "../modules/settings.h"

namespace physicell {

struct cell_definition;
struct environment;
class cell;

class Hypothesis_Rule
{
private:
	std::unordered_map<std::string, biofvm::index_t> signals_map;

public:
	std::string cell_type;
	cell_definition* pCell_Definition;

	std::string behavior;
	biofvm::real_t base_value;
	biofvm::real_t max_value;
	biofvm::real_t min_value;

	std::vector<std::string> signals;
	std::vector<bool> responses;
	std::vector<biofvm::real_t> half_maxes;
	std::vector<biofvm::real_t> hill_powers;
	std::vector<bool> applies_to_dead_cells;

	std::vector<std::string> up_signals;
	std::vector<biofvm::real_t> up_half_maxes;
	std::vector<biofvm::real_t> up_hill_powers;
	std::vector<bool> up_applies_to_dead_cells;

	std::vector<std::string> down_signals;
	std::vector<biofvm::real_t> down_half_maxes;
	std::vector<biofvm::real_t> down_hill_powers;
	std::vector<bool> down_applies_to_dead_cells;

	Hypothesis_Rule(); // done

	void sync_to_cell_definition(cell_definition* pCD, environment& e);	 // done
	void sync_to_cell_definition(std::string cell_name, environment& e); // done

	void add_signal(std::string signal, biofvm::real_t half_max, biofvm::real_t hill_power,
					std::string response);					   // done
	void add_signal(std::string signal, std::string response); // done

	biofvm::real_t evaluate(std::vector<biofvm::real_t> signal_values, bool dead); // done
	biofvm::real_t evaluate(std::vector<biofvm::real_t> signal_values);			   // done
	biofvm::real_t evaluate(cell* pCell, environment& e);						   // done
	void apply(cell* pCell, environment& e);									   // done

	biofvm::index_t find_signal(std::string name); // done

	void set_half_max(std::string, biofvm::real_t hm);	  // done
	void set_hill_power(std::string, biofvm::real_t hp);  // done
	void set_response(std::string, std::string response); // done

	void reduced_display(std::ostream& os);	 // done
	void display(std::ostream& os);			 // done
	void detailed_display(std::ostream& os); // done

	void English_display(std::ostream& os);
	void English_display_HTML(std::ostream& os);
	void English_detailed_display(std::ostream& os);
	void English_detailed_display_HTML(std::ostream& os);
};

class Hypothesis_Ruleset
{
private:
	std::unordered_map<std::string, Hypothesis_Rule*> rules_map;

public:
	std::string cell_type;
	cell_definition* pCell_Definition;

	std::vector<Hypothesis_Rule> rules;

	Hypothesis_Ruleset(); // done

	Hypothesis_Rule* add_behavior(std::string behavior, biofvm::real_t min_behavior, biofvm::real_t max_behavior,
								  environment& e);						 // done
	Hypothesis_Rule* add_behavior(std::string behavior, environment& e); // done

	// ease of access functions

	Hypothesis_Rule* find_behavior(std::string name); // done
	Hypothesis_Rule& operator[](std::string name);	  // done

	void apply(cell* pCell, environment& e);

	void sync_to_cell_definition(cell_definition* pCD, environment& e);	 // done
	void sync_to_cell_definition(std::string cell_name, environment& e); // done

	void display(std::ostream& os);			 // done
	void detailed_display(std::ostream& os); // done
};

// access

Hypothesis_Ruleset& access_ruleset(cell_definition* pCD);
Hypothesis_Ruleset* find_ruleset(cell_definition* pCD);

// initializing for all cell definitions

void intialize_hypothesis_rulesets(environment& e);

// adding and editing rules (easy eccess)

void add_rule(std::string cell_type, std::string signal, std::string behavior, std::string response, environment& e);
void add_rule(std::string cell_type, std::string signal, std::string behavior, std::string response, bool use_for_dead,
			  environment& e);

void set_hypothesis_parameters(std::string cell_type, std::string signal, std::string behavior, biofvm::real_t half_max,
							   biofvm::real_t hill_power, environment& e);
void set_behavior_parameters(std::string cell_type, std::string behavior, biofvm::real_t min_value,
							 biofvm::real_t max_value, environment& e);
void set_behavior_parameters(std::string cell_type, std::string behavior, biofvm::real_t min_value,
							 biofvm::real_t base_value, biofvm::real_t max_value, environment& e);

void set_behavior_base_value(std::string cell_type, std::string behavior, biofvm::real_t base_value, environment& e);
void set_behavior_min_value(std::string cell_type, std::string behavior, biofvm::real_t min_value, environment& e);
void set_behavior_max_value(std::string cell_type, std::string behavior, biofvm::real_t max_value, environment& e);

// display

void display_hypothesis_rulesets(std::ostream& os, environment& e);
void detailed_display_hypothesis_rulesets(std::ostream& os, environment& e);

// applying to a cell

void apply_ruleset(cell* pCell, environment& e);
void rule_phenotype_function(cell* pCell, biofvm::real_t dt, environment& e);


// parsing to / from CSV

// split a line of a CSV into a vector of strings based on the delimiter

void split_csv(std::string input, std::vector<std::string>& output, char delim);
void spaces_to_underscore(std::string& str);
std::string convert_bool_to_response(bool input);

void parse_csv_rule_v0(std::vector<std::string> input, environment& e); // parse a tokenized string (vector of strings)
void parse_csv_rule_v0(std::string input, environment& e);	   // parse a single string (a single line from CSV)
void parse_csv_rules_v0(std::string filename, environment& e); // parse all rules in a CSV file

void parse_csv_rule_v1(std::vector<std::string> input, environment& e); // parse a tokenized string (vector of strings)
void parse_csv_rule_v1(std::string input, environment& e);	   // parse a single string (a single line from CSV)
void parse_csv_rules_v1(std::string filename, environment& e); // parse all rules in a CSV file

// experimental -- removes need for base value
void parse_csv_rule_v3(std::vector<std::string> input, environment& e); // parse a tokenized string (vector of strings)
void parse_csv_rule_v3(std::string input, environment& e);	   // parse a single string (a single line from CSV)
void parse_csv_rules_v3(std::string filename, environment& e); // parse all rules in a CSV file

void parse_rules_from_pugixml(PhysiCell_Settings& settings, const pugi::xml_node& config_root, environment& e);

// needs fixing March 2023 // probably deprecate
void parse_rules_from_parameters_v0(environment& e);

std::string csv_strings_to_English(std::vector<std::string> strings, bool include_cell_header);
std::string csv_strings_to_English_v1(std::vector<std::string> strings, bool include_cell_header);
std::string csv_strings_to_English_v3(std::vector<std::string> strings, bool include_cell_header);

std::string csv_strings_to_English_HTML(std::vector<std::string> strings, bool include_cell_header);

// v1, v2, and v0?
void export_rules_csv_v0(std::string filename, environment& e);
void export_rules_csv_v1(std::string filename, environment& e);
void export_rules_csv_v3(std::string filename, environment& e);

// streamed outputs in human-readable format

void stream_annotated_English_rules(std::ostream& os, environment& e);
void stream_annotated_detailed_English_rules(std::ostream& os, environment& e);
void save_annotated_English_rules(PhysiCell_Settings& settings, environment& e);
void save_annotated_English_rules_HTML(PhysiCell_Settings& settings, environment& e);

void stream_annotated_English_rules_HTML(std::ostream& os, environment& e);
void stream_annotated_detailed_English_rules_HTML(std::ostream& os, environment& e);
void save_annotated_detailed_English_rules(PhysiCell_Settings& settings, environment& e);
void save_annotated_detailed_English_rules_HTML(PhysiCell_Settings& settings, environment& e);

// add these to PhysiCell_utilities.cpp

// std::vector<biofvm::real_t> UniformInUnitDisc();
// std::vector<biofvm::real_t> UniformInUnitSphere();

// std::vector<biofvm::real_t> UniformInAnnulus(biofvm::real_t r1, biofvm::real_t r2);
// std::vector<biofvm::real_t> UniformInShell(biofvm::real_t r1, biofvm::real_t r2);

// add this to cell behaviors

// biofvm::real_t get_single_base_behavior( Cell_Definition* pCD , std::string name );


// add these to basic signaling

/*
biofvm::real_t multivariate_Hill_response_function( std::vector<biofvm::real_t> signals, std::vector<biofvm::real_t>
half_maxes , std::vector<biofvm::real_t> hill_powers );

biofvm::real_t multivariate_linear_response_function( std::vector<biofvm::real_t> signals, std::vector<biofvm::real_t>
min_thresholds , std::vector<biofvm::real_t> max_thresholds );

std::vector<biofvm::real_t> linear_response_to_Hill_parameters( biofvm::real_t s0, biofvm::real_t s1 );
std::vector<biofvm::real_t> Hill_response_to_linear_parameters( biofvm::real_t half_max , biofvm::real_t Hill_power );
*/

void setup_cell_rules(PhysiCell_Settings& settings, const pugi::xml_node& config_root, environment& e);

}; // namespace physicell
