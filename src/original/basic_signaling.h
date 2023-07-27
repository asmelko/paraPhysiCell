#pragma once

#include <vector>


namespace physicell {

// std::vector<std::string>

double Hill_response_function(double s, double half_max, double hill_power); // done
// increases from 0 (at s_min) to 1 (at s_max)
double linear_response_function(double s, double s_min, double s_max); // done
// decreases from 1 (at s_min) to 0 (at s_max)
double decreasing_linear_response_function(double s, double s_min, double s_max); // done


double multivariate_Hill_response_function(std::vector<double> signals, std::vector<double> half_maxes,
										   std::vector<double> hill_powers);
double multivariate_linear_response_function(std::vector<double> signals, std::vector<double> min_thresholds,
											 std::vector<double> max_thresholds);

std::vector<double> linear_response_to_Hill_parameters(double s0, double s1);
std::vector<double> Hill_response_to_linear_parameters(double half_max, double Hill_power);


double interpolate_behavior(double base_value, double max_changed_value, double response);

// signal increases/decreases parameter
// options: hill power
// options: half max

class Integrated_Signal
{
private:
public:
	double base_activity;
	double max_activity;

	std::vector<double> promoters;
	std::vector<double> promoter_weights;
	double promoters_Hill;
	double promoters_half_max;

	std::vector<double> inhibitors;
	std::vector<double> inhibitor_weights;
	double inhibitors_Hill;
	double inhibitors_half_max;

	Integrated_Signal();
	void reset(void);

	void add_signal(char signal_type, double signal, double weight);
	void add_signal(char signal_type, double signal);

	double compute_signal(void);
};


}; // namespace physicell
