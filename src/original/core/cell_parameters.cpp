#include "cell_parameters.h"

#include "constants.h"

using namespace physicell;

Cell_Parameters::Cell_Parameters()
{
	o2_hypoxic_threshold = 15.0; // HIF-1alpha at half-max around 1.5-2%, and tumors often are below 2%
	o2_hypoxic_response = 8.0;	 // genomic / proteomic changes observed at 7-8 mmHg
	o2_hypoxic_saturation = 4.0; // maximum HIF-1alpha at 0.5% o2 (McKeown)

	o2_necrosis_threshold = 5.0;
	o2_necrosis_max = 2.5;

	o2_proliferation_threshold = 5.0;	 // assume no proliferation at same level as starting necrosis
	o2_proliferation_saturation = 160.0; // 5% = 38, 21% = 160 mmHg
	o2_reference = 160.0;				 // assume all was measured in normoxic 21% o2

	pReference_live_phenotype = NULL; // reference live (usually physioxic) phenotype

	// necrosis parameters

	max_necrosis_rate = 1.0 / (6.0 * 60.0); // assume cells survive 6 hours in very low oxygen
	necrosis_type = constants::deterministic_necrosis;
}
