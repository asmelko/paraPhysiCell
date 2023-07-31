#pragma once

#include <BioFVM/types.h>

namespace physicell {
namespace constants {

constexpr biofvm::index_t deterministic_necrosis = 0;
constexpr biofvm::index_t stochastic_necrosis = 1;

// currently recognized cell cycle models
constexpr biofvm::index_t advanced_Ki67_cycle_model = 0;
constexpr biofvm::index_t basic_Ki67_cycle_model = 1;
constexpr biofvm::index_t flow_cytometry_cycle_model = 2;
constexpr biofvm::index_t live_apoptotic_cycle_model = 3;
constexpr biofvm::index_t total_cells_cycle_model = 4;
constexpr biofvm::index_t live_cells_cycle_model = 5;
constexpr biofvm::index_t flow_cytometry_separated_cycle_model = 6;
constexpr biofvm::index_t cycling_quiescent_model = 7;

// currently recognized death models
constexpr biofvm::index_t apoptosis_death_model = 100;
constexpr biofvm::index_t necrosis_death_model = 101;
constexpr biofvm::index_t autophagy_death_model = 102;

constexpr biofvm::index_t custom_cycle_model = 9999;

// currently recognized cell cycle and death phases
// cycle phases
constexpr biofvm::index_t Ki67_positive_premitotic = 0;
constexpr biofvm::index_t Ki67_positive_postmitotic = 1;
constexpr biofvm::index_t Ki67_positive = 2;
constexpr biofvm::index_t Ki67_negative = 3;
constexpr biofvm::index_t G0G1_phase = 4;
constexpr biofvm::index_t G0_phase = 5;
constexpr biofvm::index_t G1_phase = 6;
constexpr biofvm::index_t G1a_phase = 7;
constexpr biofvm::index_t G1b_phase = 8;
constexpr biofvm::index_t G1c_phase = 9;
constexpr biofvm::index_t S_phase = 10;
constexpr biofvm::index_t G2M_phase = 11;
constexpr biofvm::index_t G2_phase = 12;
constexpr biofvm::index_t M_phase = 13;
constexpr biofvm::index_t live = 14;

constexpr biofvm::index_t G1pm_phase = 15;
constexpr biofvm::index_t G1ps_phase = 16;

constexpr biofvm::index_t cycling = 17;
constexpr biofvm::index_t quiescent = 18;

constexpr biofvm::index_t custom_phase = 9999;

// death phases
constexpr biofvm::index_t apoptotic = 100;
constexpr biofvm::index_t necrotic_swelling = 101;
constexpr biofvm::index_t necrotic_lysed = 102;
constexpr biofvm::index_t necrotic = 103;
constexpr biofvm::index_t debris = 104;

} // namespace constants
} // namespace physicell
