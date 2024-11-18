#include "various_output.h"

#include <cmath>
#include <cstring>

#include "../../environment.h"
#include "../core/constants.h"
#include "settings.h"
#include "timer.h"

namespace physicell {

void display_simulation_status(std::ostream& os, environment& e, PhysiCell_Settings& settings)
{
	os << "current simulated time: " << e.current_time << " " << settings.time_units << " (max: " << settings.max_time
	   << " " << settings.time_units << ")" << std::endl;

	os << "total agents: " << e.m.agents->agents_count() << std::endl;

	os << "interval wall time: ";
	TOC();
	display_stopwatch_value(os, stopwatch_value());
	os << std::endl;
	TIC();

	os << "total wall time: ";
	RUNTIME_TOC();
	display_stopwatch_value(os, runtime_stopwatch_value());
	os << std::endl << std::endl;

	return;
}

int writePov(environment& e, PhysiCell_Settings& settings, double timepoint, double scale)
{
	static int TUMOR_TYPE = 0;
	static int VESSEL_TYPE = 1;

	std::string filename;
	filename.resize(1024);
	//	sprintf( (char*) filename.c_str() , "output//cells_%i.pov" , (int)round(timepoint) );
	sprintf((char*)filename.c_str(), "%s/cells_%i.pov", settings.folder.c_str(), (int)round(timepoint));
	std::ofstream povFile(filename.c_str(), std::ofstream::out);
	povFile << "#include \"colors.inc\" \n";
	povFile << "#include \"header.inc\" \n";

	for (int i = 0; i < e.m.agents->agents_count(); i++)
	{
		std::string _nameCore;

		if (e.get_container().get_at(i)->phenotype.cycle.pCycle_Model)
		{
			int code = e.get_container().get_at(i)->phenotype.cycle.current_phase().code;
			if (code == constants::Ki67_positive_premitotic || code == constants::Ki67_positive_postmitotic
				|| code == constants::Ki67_positive || code == constants::Ki67_negative || code == constants::live)
				_nameCore = "LIVE";
			else if (code == constants::apoptotic)
				_nameCore = "APOP";
			else if (code == constants::necrotic_swelling || code == constants::necrotic_lysed
					 || code == constants::necrotic)
				_nameCore = "NEC";
			else if (code == constants::debris)
				_nameCore = "DEBR";
			else
				_nameCore = "MISC";
		}
		else if (e.get_container().get_at(i)->type == TUMOR_TYPE)
			_nameCore = "LIVE";
		else if (e.get_container().get_at(i)->type == VESSEL_TYPE)
			_nameCore = "ENDO";
		else
			_nameCore = "MISC";
		std::string center = std::string("<") + std::to_string(e.get_container().get_at(i)->get_position()[0] / scale)
							 + "," + std::to_string(e.get_container().get_at(i)->get_position()[1] / scale) + ","
							 + std::to_string(e.get_container().get_at(i)->get_position()[2] / scale) + ">";
		std::string core = std::string("sphere {\n\t") + center + "\n\t "
						   + std::to_string(e.get_container().get_at(i)->phenotype.geometry.radius() / scale)
						   + "\n\t FinishMacro ( " + center + "," + _nameCore + "Finish," + _nameCore + "*1)\n}\n";
		povFile << core;
	}

	povFile << "#include \"footer.inc\" \n";
	povFile.close();
	return 0;
}

int writeCellReport(environment& e, PhysiCell_Settings& settings, double timepoint)
{
	std::string filename;
	filename.resize(1024);
	//	sprintf( (char*) filename.c_str() , "output//cells_%i.txt" , (int)round(timepoint) );
	sprintf((char*)filename.c_str(), "%s/cells_%i.txt", settings.folder.c_str(), (int)round(timepoint));
	std::ofstream povFile(filename.c_str(), std::ofstream::out);
	povFile << "\tID\tx\ty\tz\tradius\tvolume_total\tvolume_nuclear_fluid\tvolume_nuclear_solid\tvolume_cytoplasmic_"
			   "fluid\tvolume_cytoplasmic_solid\tvolume_calcified_fraction\tphenotype\telapsed_time\n";
	int phenotype_code;
	for (int i = 0; i < e.m.agents->agents_count(); i++)
	{
		auto cell = e.get_container().get_at(i);
		phenotype_code = cell->phenotype.cycle.current_phase().code;
		// phenotype_code =
		// phases.size()>0?cell->phenotype.cycle.phases[cell->phenotype.current_phase_index].code:-1;
		povFile << i << "\t" << cell->id << "\t" << cell->get_position()[0] << "\t" << cell->get_position()[1] << "\t"
				<< cell->get_position()[2] << "\t";
		povFile << cell->phenotype.geometry.radius() << "\t" << cell->phenotype.volume.total() << "\t"
				<< cell->phenotype.volume.nuclear_fluid() << "\t" << cell->phenotype.volume.nuclear_solid() << "\t"
				<< cell->phenotype.volume.cytoplasmic_fluid() << "\t" << cell->phenotype.volume.cytoplasmic_solid()
				<< "\t" << cell->phenotype.volume.calcified_fraction() << "\t" << phenotype_code <<
			// "\t"<< cell->phenotype.cycle.phases[cell->phenotype.current_phase_index].elapsed_time
			// <<std::endl;
			"\t" << cell->phenotype.cycle.data.elapsed_time_in_phase << std::endl;
	}
	povFile.close();
	return 0;
}

void log_output(double t, int output_index, environment& e, PhysiCell_Settings& settings, std::ofstream& report_file)
{
	double scale = 1000;
	int num_new_cells = 0;
	int num_deaths = 0;
	//	std::cout << "current simulated time: " << t   << " minutes " << std::endl;
	//	std::cout << "interval wall time: ";
	//	BioFVM::TOC();
	//	BioFVM::display_stopwatch_value( std::cout , BioFVM::stopwatch_value() );
	//	std::cout << std::endl;
	//	std::cout << "total wall time: ";
	//	BioFVM::RUNTIME_TOC();
	//	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() );
	//	std::cout << std::endl;

	std::cout << "time: " << t << std::endl;
	num_new_cells = t == 0 ? e.m.agents->agents_count() : e.divisions_count;
	num_deaths = e.deaths_count;
	std::cout << "total number of agents (newly born, deaths): " << e.m.agents->agents_count() << "(" << num_new_cells
			  << ", " << num_deaths << ")" << std::endl;
	report_file << t << "\t" << e.m.agents->agents_count() << "\t" << num_new_cells << "\t" << num_deaths << "\t"
				<< stopwatch_value() << std::endl;
	//	BioFVM::TIC();

	e.divisions_count = 0;
	e.deaths_count = 0;
	writePov(e, settings, t, scale);
	writeCellReport(e, settings, t);
	std::string filename;
	filename.resize(1024, '\0');
	sprintf((char*)filename.c_str(), "output%08d.mat", output_index);
	filename.resize(strlen(filename.c_str()));
	// std::cout << "\tWriting to file " << filename << " ... " << std::endl;
	// microenvironment.write_to_matlab( filename );

	return;
}

} // namespace physicell
