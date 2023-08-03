#pragma once

#include <pugixml.hpp>
#include <string>

namespace biofvm {
struct microenvironment;
}

namespace physicell {

struct environment;

// void add_PhysiCell_cell_to_open_xml_pugi(  pugi::xml_document& xml_dom, Cell& C ); // not implemented -- future
// edition
void add_PhysiCell_cells_to_open_xml_pugi(pugi::xml_document& xml_dom, std::string filename_base,
										  biofvm::microenvironment& M);
void add_PhysiCell_to_open_xml_pugi(pugi::xml_document& xml_dom, std::string filename_base,
									double current_simulation_time, biofvm::microenvironment& M);


void save_PhysiCell_to_MultiCellDS_xml_pugi(std::string filename_base, environment& e);


/* V2 functions */

/*
void add_PhysiCell_cell_to_open_xml_pugi_v2(  pugi::xml_document& xml_dom, Cell& C ); // not implemented -- future
edition void add_PhysiCell_cells_to_open_xml_pugi_v2( pugi::xml_document& xml_dom, std::string filename_base,
biofvm::microenvironment& M  ); void add_PhysiCell_to_open_xml_pugi_v2( pugi::xml_document& xml_dom , std::string
filename_base, double current_simulation_time , biofvm::microenvironment& M );

void save_PhysiCell_to_MultiCellDS_xml_pugi_v2( std::string filename_base , biofvm::microenvironment& M , double
current_simulation_time);
*/

void add_PhysiCell_cells_to_open_xml_pugi_v2(pugi::xml_document& xml_dom, std::string filename_base, environment& e);
void save_PhysiCell_to_MultiCellDS_v2(std::string filename_base, environment& e);
void write_neighbor_graph(std::string filename, environment& e);
void write_attached_cells_graph(std::string filename, environment& e);

}; // namespace physicell
