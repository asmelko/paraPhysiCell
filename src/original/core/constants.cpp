#include "constants.h"

#include <string>
#include <unordered_map>

using namespace biofvm;

namespace physicell {
std::unordered_map<const char*, index_t> cycle_model_codes = {
	{ "Ki67 (advanced)", constants::advanced_Ki67_cycle_model },
	{ "Ki67 (basic)", constants::basic_Ki67_cycle_model },
	{ "Flow cytometry model (basic)", constants::flow_cytometry_cycle_model },
	// { ,constants::live_apoptotic_cycle_model}, // not implemented
	// { ,constants::total_cells_cycle_model}, // not implemented
	{ "Live", constants::live_cells_cycle_model },
	{ "Flow cytometry model (separated)", constants::flow_cytometry_separated_cycle_model },
	{ "Cycling-Quiescent model", constants::cycling_quiescent_model },

	// currently recognized death models
	{ "Apoptosis", constants::apoptosis_death_model },
	{ "Necrosis", constants::necrosis_death_model },
	// { ,constants::autophagy_death_model}, // not implemented

	{ "ki67 (advanced)", constants::advanced_Ki67_cycle_model },
	{ "ki67 (basic)", constants::basic_Ki67_cycle_model },
	{ "flow cytometry model (basic)", constants::flow_cytometry_cycle_model },
	{ "live", constants::live_cells_cycle_model },
	{ "flow cytometry model (separated)", constants::flow_cytometry_separated_cycle_model },
	{ "cycling-quiescent model", constants::cycling_quiescent_model },
	{ "apoptosis", constants::apoptosis_death_model },
	{ "necrosis", constants::necrosis_death_model }

};

index_t find_cycle_model_code(const std::string& model_name)
{
	auto search = cycle_model_codes.find(model_name.c_str());
	if (search == cycle_model_codes.end())
	{
		return -1;
	}
	else
	{
		return search->second;
	}
}
} // namespace physicell
