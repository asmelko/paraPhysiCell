#pragma once

#include <pugixml.hpp>
#include <string>

namespace biofvm {

extern std::string MultiCellDS_version_string;

extern std::string MultiCellDS_clinical_snapshot_type_string;
extern int MultiCellDS_clinical_snapshot_code;

extern std::string MultiCellDS_experimental_snapshot_type_string;
extern int MultiCellDS_experimental_snapshot_code;

extern std::string MultiCellDS_simulation_snapshot_type_string;
extern int MultiCellDS_simulation_snapshot_code;

extern std::string MultiCellDS_digital_cell_line_type_string;
extern int MultiCellDS_digital_cell_line_code;

/* options */

extern bool save_mesh_as_matlab;
extern bool save_density_data_as_matlab;
extern bool save_cells_as_custom_matlab;
extern bool save_cell_data;


struct microenvironment;

extern pugi::xml_document biofvm_doc;

class Person_Metadata
{
private:
	bool is_empty;

public:
	std::string type; // author, creator, user, curator
	std::string surname;
	std::string given_names;
	std::string email;
	std::string URL;
	std::string organization;
	std::string department;
	std::string ORCID;

	Person_Metadata();
	void display_information(std::ostream& os);
	void insert_in_open_xml_pugi(pugi::xml_node& insert_here);
};

class Citation_Metadata
{
private:
public:
	std::string DOI;
	std::string PMID;
	std::string PMCID;
	std::string text;
	std::string notes;
	std::string URL;

	Citation_Metadata();
	void display_information(std::ostream& os);
	void insert_in_open_xml_pugi(pugi::xml_node& insert_here);
};

class Software_Metadata
{
private:
public:
	// basic program information
	std::string program_name;
	std::string program_version;
	std::string program_URL;

	Person_Metadata creator;
	Person_Metadata user;
	Citation_Metadata citation;

	Software_Metadata();

	void display_information(std::ostream& os);
	void insert_in_open_xml_pugi(pugi::xml_node& insert_here);
};

class MultiCellDS_Metadata
{
private:
public:
	std::string MultiCellDS_type;

	Software_Metadata program;
	Citation_Metadata data_citation;

	// scientific information
	std::string spatial_units;
	std::string time_units;
	std::string runtime_units;
	double current_time;
	double current_runtime;

	std::string description; // any optional text -- not implemented

	MultiCellDS_Metadata();
	void display_information(std::ostream& os);
	void sync_to_microenvironment(microenvironment& M);
	void restart_runtime(void);

	void add_to_open_xml_pugi(double current_simulation_time, pugi::xml_document& xml_dom);
};

extern MultiCellDS_Metadata BioFVM_metadata;

/* setting up the main MultiCellDS tree structure */

void add_MultiCellDS_main_structure_to_open_xml_pugi(pugi::xml_document& xml_dom);

/* set options */

void set_save_biofvm_mesh_as_matlab(bool newvalue);				// default: true
void set_save_biofvm_data_as_matlab(bool newvalue);				// default: true
void set_save_biofvm_cell_data(bool newvalue);					// default: true
void set_save_biofvm_cell_data_as_custom_matlab(bool newvalue); // default: true

/* writing parts of BioFVM to a MultiCellDS file */

void add_BioFVM_substrates_to_open_xml_pugi(pugi::xml_document& xml_dom, std::string filename_base,
											microenvironment& M);
// void add_BioFVM_basic_agent_to_open_xml_pugi(  pugi::xml_document& xml_dom, Basic_Agent& BA ); // not implemented --
// future edition
void add_BioFVM_agents_to_open_xml_pugi(pugi::xml_document& xml_dom, std::string filename_base, microenvironment& M);
void add_BioFVM_to_open_xml_pugi(pugi::xml_document& xml_dom, std::string filename_base, double current_simulation_time,
								 microenvironment& M);

void save_BioFVM_to_MultiCellDS_xml_pugi(std::string filename_base, microenvironment& M,
										 double current_simulation_time);

// /* beta in PhysiCell 1.11.0 */

// bool read_microenvironment_from_matlab(std::string mat_filename);

// /* future / not yet supported */

// void read_BioFVM_from_open_xml_pugi(pugi::xml_document& xml_dom, std::string filename_base,
// 									double& current_simulation_time, microenvironment& M);
// void read_BioFVM_to_MultiCellDS_xml_pugi(std::string filename_base, microenvironment& M,
// 										 double& current_simulation_time);

// /* partly-implemented code snippets -- not to be used as of March 2016 */

// // functions to read multiscale_microenvironment from MultiCellDS file (requires pugixml)
// void read_microenvironment_from_MultiCellDS_xml(microenvironment& M_destination, std::string filename);
// void read_microenvironment_from_MultiCellDS_xml(microenvironment& M_destination, pugi::xml_document& xml_dom);

}; // namespace biofvm
