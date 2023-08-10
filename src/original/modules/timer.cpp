#include "timer.h"

#include <chrono>
#include <cmath>

namespace physicell {

std::chrono::steady_clock::time_point tic_time;
std::chrono::steady_clock::time_point toc_time;
std::chrono::steady_clock::time_point program_tic_time;
std::chrono::steady_clock::time_point program_toc_time;

double total_tictoc_time = 0.0;

void TIC(void) { tic_time = std::chrono::steady_clock::now(); }

void TOC(void)
{
	toc_time = std::chrono::steady_clock::now();
	total_tictoc_time += stopwatch_value();
}

void RUNTIME_TIC(void) { program_tic_time = std::chrono::steady_clock::now(); }

void RUNTIME_TOC(void) { program_toc_time = std::chrono::steady_clock::now(); }

double stopwatch_value(void)
{
	static std::chrono::duration<double> time_span;
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(toc_time - tic_time);
	return time_span.count();
}

double runtime_stopwatch_value(void)
{
	static std::chrono::duration<double> time_span;
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(program_toc_time - program_tic_time);
	return time_span.count();
}

std::string format_stopwatch_value(double dIn)
{
	std::string output;
	output.resize(1024);
	int nDays = (int)floor((double)(dIn / (60.0 * 60.0 * 24.0)));
	int nHours = (int)floor((double)((dIn - nDays * 60 * 60 * 24) / (60.0 * 60.0)));
	int nMinutes = (int)floor((double)((dIn - nDays * 60 * 60 * 24 - nHours * 60 * 60) / (60.0)));
	double dSeconds = dIn - nDays * 60.0 * 60.0 * 24.0 - nHours * 60.0 * 60.0 - nMinutes * 60.0;

	sprintf((char*)output.c_str(), "%d days, %d hours, %d minutes, and %2.4f seconds", nDays, nHours, nMinutes,
			dSeconds);

	return output;
}

void display_stopwatch_value(std::ostream& os, double dIn)
{
	int nDays = (int)floor((double)(dIn / (60.0 * 60.0 * 24.0)));
	int nHours = (int)floor((double)((dIn - nDays * 60 * 60 * 24) / (60.0 * 60.0)));
	int nMinutes = (int)floor((double)((dIn - nDays * 60 * 60 * 24 - nHours * 60 * 60) / (60.0)));
	double dSeconds = dIn - nDays * 60.0 * 60.0 * 24.0 - nHours * 60.0 * 60.0 - nMinutes * 60.0;

	os << nDays << " days, " << nHours << " hours, " << nMinutes << " minutes, and " << dSeconds << " seconds ";
	return;
}

} // namespace physicell