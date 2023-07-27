#include "basic_signaling.h"

#include <cmath>

namespace physicell {

Integrated_Signal::Integrated_Signal()
{
	base_activity = 0.0;
	max_activity = 1.0;

	promoters.clear();
	promoter_weights.clear();

	promoters_half_max = 0.1;
	promoters_Hill = 4;

	inhibitors.clear();
	inhibitor_weights.clear();

	inhibitors_half_max = 0.1;
	inhibitors_Hill = 4;

	return;
}

void Integrated_Signal::reset(void)
{
	promoters.clear();
	promoter_weights.clear();

	inhibitors.clear();
	inhibitor_weights.clear();
	return;
}

double Integrated_Signal::compute_signal(void)
{
	double pr = 0.0;
	double w = 0.0;
	for (std::size_t k = 0; k < promoters.size(); k++)
	{
		pr += promoters[k];
		w += promoter_weights[k];
	}
	w += 1e-16;
	pr /= w;

	double inhib = 0.0;
	w = 0.0;
	for (std::size_t k = 0; k < inhibitors.size(); k++)
	{
		inhib += inhibitors[k];
		w += inhibitor_weights[k];
	}
	w += 1e-16;
	inhib /= w;

	double Pn = pow(pr, promoters_Hill);
	double Phalf = pow(promoters_half_max, promoters_Hill);

	double In = pow(inhib, inhibitors_Hill);
	double Ihalf = pow(inhibitors_half_max, inhibitors_Hill);

	double P = Pn / (Pn + Phalf);
	double I = 1.0 / (In + Ihalf);

	double output = max_activity;
	output -= base_activity; //(max-base)
	output *= P;			 // (max-base)*P
	output += base_activity; // base + (max-base)*P
	output *= I;			 // (base + (max-base)*P)*I;

	return output;
};

void Integrated_Signal::add_signal(char signal_type, double signal, double weight)
{
	if (signal_type == 'P' || signal_type == 'p')
	{
		promoters.push_back(signal);
		promoter_weights.push_back(weight);
		return;
	}
	if (signal_type == 'I' || signal_type == 'i')
	{
		inhibitors.push_back(signal);
		inhibitor_weights.push_back(weight);
		return;
	}
	return;
}

void Integrated_Signal::add_signal(char signal_type, double signal) { return add_signal(signal_type, signal, 1.0); }

double Hill_response_function(double s, double half_max, double hill_power)
{
	// newer. only one expensive a^b operation. 45% less computationl expense.

	// give an early exit possibility to cut cost on "empty" rules
	if (s < 1e-16) // maybe also try a dynamic threshold: 0.0001 * half_max
	{
		return 0.0;
	}

	// operations to reduce a^b operations and minimize hidden memory allocation / deallocation / copy operations.
	// Hill = (s/half_max)^hill_power / ( 1 + (s/half_max)^hill_power  )
	double temp = s;					  // s
	temp /= half_max;					  // s/half_max
	double temp1 = pow(temp, hill_power); // (s/half_max)^h
	temp = temp1;						  // (s/half_max)^h
	temp += 1;							  // (1+(s/half_max)^h );
	temp1 /= temp;						  // (s/half_max)^h / ( 1 + s/half_max)^h)
	return temp1;
}

double linear_response_function(double s, double s_min, double s_max)
{
	if (s <= s_min)
	{
		return 0.0;
	}
	if (s >= s_max)
	{
		return 1.0;
	}
	s -= s_min;		// overwrite s with s - s_min
	s_max -= s_min; // overwrite s_max with s_max - s_min
	s /= s_max;		// now we have (s-s_min)/(s_max-s_min
	return s;
}

double decreasing_linear_response_function(double s, double s_min, double s_max)
{
	if (s <= s_min)
	{
		return 1.0;
	}
	if (s >= s_max)
	{
		return 0.0;
	}
	// (smax-s)/(smax-smin);
	// = -(s-smax)/(smax-smin)
	s -= s_max;		// replace s by s-s_max
	s_max -= s_min; // replace s_max = s_max - s_min
	s /= s_max;		// this is (s-s_max)/(s_max-s_min)
	s *= -1;		// this is (s_max-s)/(s_max-s_min)
	return s;
}

double interpolate_behavior(double base_value, double max_changed_value, double response)
{
	double output = max_changed_value; // bM
	output -= base_value;			   // (bM-b0);
	output *= response;				   // R*(bM-b0);
	output += base_value;			   // b0 + (bM-b0)*R;
	return output;
}

double multivariate_Hill_response_function(std::vector<double> signals, std::vector<double> half_maxes,
										   std::vector<double> hill_powers)
{
	double temp1 = 0.0;
	double temp2 = 0.0;
	double temp3 = 0.0;
	// create the generalized (s^h), stored in temp1;
	for (std::size_t j = 0; j < signals.size(); j++)
	{
		temp2 = signals[j];					// s
		temp2 /= half_maxes[j];				// s/s_half
		temp3 = pow(temp2, hill_powers[j]); // (s/s_half)^h
		temp1 += temp3;
	}
	temp2 = temp1;	// numerator (S^h)
	temp1 += 1.0;	// denominator (1+S^h)
	temp2 /= temp1; // numerator/denominator = S^h / (1+S^h)
	return temp2;
}

double multivariate_linear_response_function(std::vector<double> signals, std::vector<double> min_thresholds,
											 std::vector<double> max_thresholds)
{
	double output = 0.0;

	for (std::size_t j = 0; j < signals.size(); j++)
	{
		output += linear_response_function(signals[j], min_thresholds[j], max_thresholds[j]);
	}

	if (output > 1.0)
	{
		return 1.0;
	}

	return output;
}

std::vector<double> linear_response_to_Hill_parameters(double s0, double s1)
{
	static double tol = 0.1;
	static double param1 = (1 - tol) / tol;
	static double param2 = log(param1);

	// half max, then hill power
	double hm = 0.5 * (s0 + s1);

	// hp so that H(s1) ~ (1-tol)
	double hp = round(param2 / log(s1 / hm));

	std::vector<double> output = { hm, hp };

	return output;
}

std::vector<double> Hill_response_to_linear_parameters(double half_max, double Hill_power)
{
	static double tol = 0.1;
	static double param1 = (1 - tol) / tol;
	double param2 = pow(param1, 1.0 / Hill_power);

	// s1 such that H(s1) ~ (1-tol)
	double s1 = half_max * param2;

	// s0 for symmetry
	double s0 = 2 * half_max - s1;
	if (s0 < 0)
	{
		s0 = 0.0;
	}

	std::vector<double> output = { s0, s1 };

	return output;
}



}; // namespace physicell
