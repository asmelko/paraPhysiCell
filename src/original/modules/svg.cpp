#include "svg.h"

#include <cstring>

namespace physicell {

bool Write_SVG_start(std::ostream& os, double width, double height)
{
	os << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << std::endl
	   << "<!-- Created with PhysiCell (http://PhysiCell.MathCancer.org/) -->" << std::endl;

	os << "<svg " << std::endl
	   << " xmlns:dc=\"http://purl.org/dc/elements/1.1/\" " << std::endl
	   << " xmlns:cc=\"http://creativecommons.org/ns#\" " << std::endl
	   << " xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" " << std::endl
	   << " xmlns:svg=\"http://www.w3.org/2000/svg\" " << std::endl
	   << " xmlns=\"http://www.w3.org/2000/svg\" " << std::endl
	   << " version=\"1.1\" " << std::endl
	   << " width=\"" << width << "\" " << std::endl
	   << " height=\"" << height << "\" " << std::endl
	   << " id=\"svg2\">" << std::endl;

	return true;
}

bool Write_SVG_end(std::ostream& os)
{
	os << "</svg>" << std::endl;
	return true;
}

bool Write_SVG_text(std::ostream& os, const char* str, double position_x, double position_y, double font_size,
					const char* color, const char* font)
{
	os << "  <text x=\"" << position_x << "\" y=\"" << position_y << "\"" << std::endl
	   << "   font-family=\"" << font << "\" font-size=\"" << font_size << "\" fill=\"" << color << "\" >" << std::endl
	   << "   " << str << std::endl
	   << "  </text>" << std::endl;
	return true;
}

void Write_SVG_text(std::ostream& os, const char* str, double position_x, double position_y, double font_size,
					const char* color, const char* font, double rotation)
{
	double text_width = font_size * std::strlen(str) / 2.0; // estimate the width of the text
	double text_height = font_size / 2.0;					// estimate the height of the text

	double center_x = position_x + text_width / 2.0;
	double center_y = position_y + text_height / 2.0;

	os << "<text x=\"" << position_x << "\" y=\"" << position_y << "\" font-size=\"" << font_size << "\" fill=\""
	   << color << "\" font-family=\"" << font << "\" transform=\"rotate(" << rotation << " " << center_x << " "
	   << center_y << ")\">" << str << "</text>\n";
}

bool Write_SVG_circle(std::ostream& os, double center_x, double center_y, double radius, double stroke_size,
					  std::string stroke_color, std::string fill_color)
{
	os << "  <circle cx=\"" << center_x << "\" cy=\"" << center_y << "\" r=\"" << radius << "\" stroke-width=\""
	   << stroke_size << "\" stroke=\"" << stroke_color << "\" fill=\"" << fill_color << "\"/>" << std::endl;
	return true;
}


bool Write_SVG_rect(std::ostream& os, double UL_corner_x, double UL_corner_y, double width, double height,
					double stroke_size, std::string stroke_color, std::string fill_color)
{
	os << "  <rect x=\"" << UL_corner_x << "\" y=\"" << UL_corner_y << "\" width=\"" << width << "\" height=\""
	   << height << "\" stroke-width=\"" << stroke_size << "\" stroke=\"" << stroke_color << "\" fill=\"" << fill_color
	   << "\"/>" << std::endl;
	return true;
}

bool Write_SVG_line(std::ostream& os, double start_x, double start_y, double end_x, double end_y, double thickness,
					std::string stroke_color)
{
	os << "  <line x1=\"" << start_x << "\" y1=\"" << start_y << "\" x2=\"" << end_x << "\" y2=\"" << end_y << "\" "
	   << "stroke=\"" << stroke_color << "\" stroke-width=\"" << thickness << "\"/>" << std::endl;
	return true;
}

} // namespace physicell
