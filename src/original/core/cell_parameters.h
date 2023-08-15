#pragma once

#include "../../phenotype.h"

namespace physicell {

class Cell_Parameters
{
private:
public:
	// oxygen values (in mmHg) for critical phenotype changes
	biofvm::real_t o2_hypoxic_threshold;  // value at which hypoxic signaling starts
	biofvm::real_t o2_hypoxic_response;	  // value at which omics changes are observed
	biofvm::real_t o2_hypoxic_saturation; // value at which hypoxic signalign saturates
	// o2_hypoxic_saturation < o2_hypoxic_threshold

	biofvm::real_t o2_proliferation_saturation; // value at which extra o2 does not increase proliferation
	biofvm::real_t o2_proliferation_threshold;	// value at which o2 is sufficient for proliferation

	biofvm::real_t o2_reference; // physioxic reference value, in the linked reference Phenotype
	// o2_proliferation_threshold < o2_reference < o2_proliferation_saturation;

	biofvm::real_t o2_necrosis_threshold; // value at which cells start experiencing necrotic death
	biofvm::real_t o2_necrosis_max;		  // value at which necrosis reaches its maximum rate
	// o2_necrosis_max < o2_necrosis_threshold

	phenotype_t* pReference_live_phenotype; // reference live phenotype (typically physioxic)
	// Phenotype* pReference_necrotic_phenotype; // reference live phenotype (typically physioxic)

	// necrosis parameters (may evenually be moved into a reference necrotic phenotype
	biofvm::real_t max_necrosis_rate; // deprecate
	biofvm::index_t necrosis_type;	  // deprecate

	Cell_Parameters();
};

} // namespace physicell
