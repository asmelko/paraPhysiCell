#pragma once

#include <string>

namespace physicell {
FILE* write_matlab_header(unsigned int rows, unsigned int cols, std::string filename, std::string variable_name);
}
