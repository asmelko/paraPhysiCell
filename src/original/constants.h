#pragma once

#include <BioFVM/types.h>

namespace physicell {
namespace constants {

constexpr biofvm::index_t deterministic_necrosis = 0;
constexpr biofvm::index_t stochastic_necrosis = 1;

// currently recognized death models
constexpr biofvm::index_t apoptosis_death_model = 100;
constexpr biofvm::index_t necrosis_death_model = 101;
constexpr biofvm::index_t autophagy_death_model = 102;

constexpr biofvm::index_t custom_cycle_model = 9999;

// death phases
constexpr biofvm::index_t apoptotic = 100;
constexpr biofvm::index_t necrotic_swelling = 101;
constexpr biofvm::index_t necrotic_lysed = 102;
constexpr biofvm::index_t necrotic = 103;
constexpr biofvm::index_t debris = 104;

} // namespace constants
} // namespace physicell
