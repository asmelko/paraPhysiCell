#pragma once

#include <string>
#include <vector>

#include <BioFVM/types.h>

namespace physicell {
std::vector<biofvm::real_t> csv_to_vector(const std::string& value);
}
