#pragma once

#include <string>

namespace physicell {

void TIC(void);

void TOC(void);

void RUNTIME_TIC(void);

void RUNTIME_TOC(void);

double stopwatch_value(void);

double runtime_stopwatch_value(void);

std::string format_stopwatch_value(double dIn);

} // namespace physicell
