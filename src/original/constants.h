#pragma once

#include <BioFVM/types.h>

namespace physicell {
namespace constants {

// currently recognized death models
constexpr biofvm::index_t apoptosis_death_model = 100;
constexpr biofvm::index_t necrosis_death_model = 101;
constexpr biofvm::index_t autophagy_death_model = 102;

constexpr biofvm::index_t custom_cycle_model = 9999;

} // namespace constants
} // namespace physicell
