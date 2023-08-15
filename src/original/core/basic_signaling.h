#pragma once

#include <vector>

#include <BioFVM/types.h>

namespace physicell {

// std::vector<std::string>

biofvm::real_t Hill_response_function(biofvm::real_t s, biofvm::real_t half_max, biofvm::real_t hill_power); // done
// increases from 0 (at s_min) to 1 (at s_max)
biofvm::real_t linear_response_function(biofvm::real_t s, biofvm::real_t s_min, biofvm::real_t s_max); // done
// decreases from 1 (at s_min) to 0 (at s_max)
biofvm::real_t decreasing_linear_response_function(biofvm::real_t s, biofvm::real_t s_min,
												   biofvm::real_t s_max); // done


biofvm::real_t multivariate_Hill_response_function(std::vector<biofvm::real_t> signals,
												   std::vector<biofvm::real_t> half_maxes,
												   std::vector<biofvm::real_t> hill_powers);
biofvm::real_t multivariate_linear_response_function(std::vector<biofvm::real_t> signals,
													 std::vector<biofvm::real_t> min_thresholds,
													 std::vector<biofvm::real_t> max_thresholds);

std::vector<biofvm::real_t> linear_response_to_Hill_parameters(biofvm::real_t s0, biofvm::real_t s1);
std::vector<biofvm::real_t> Hill_response_to_linear_parameters(biofvm::real_t half_max, biofvm::real_t Hill_power);


biofvm::real_t interpolate_behavior(biofvm::real_t base_value, biofvm::real_t max_changed_value,
									biofvm::real_t response);

// signal increases/decreases parameter
// options: hill power
// options: half max

class Integrated_Signal
{
public:
	biofvm::real_t base_activity;
	biofvm::real_t max_activity;

	std::vector<biofvm::real_t> promoters;
	std::vector<biofvm::real_t> promoter_weights;
	biofvm::real_t promoters_Hill;
	biofvm::real_t promoters_half_max;

	std::vector<biofvm::real_t> inhibitors;
	std::vector<biofvm::real_t> inhibitor_weights;
	biofvm::real_t inhibitors_Hill;
	biofvm::real_t inhibitors_half_max;

	Integrated_Signal();
	void reset();

	void add_signal(char signal_type, biofvm::real_t signal, biofvm::real_t weight);
	void add_signal(char signal_type, biofvm::real_t signal);

	biofvm::real_t compute_signal();
};


}; // namespace physicell
