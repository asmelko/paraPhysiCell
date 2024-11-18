#pragma once

#include <iostream>
#include <string>

namespace physicell {

bool Write_SVG_start(std::ostream& os, double width, double height);
bool Write_SVG_end(std::ostream& os);

bool Write_SVG_text(std::ostream& os, const char* str, double position_x, double position_y, double font_size,
					const char* color, const char* font);
void Write_SVG_text(std::ostream& os, const char* str, double position_x, double position_y, double font_size,
					const char* color, const char* font, double rotation);

bool Write_SVG_circle(std::ostream& os, double center_x, double center_y, double radius, double stroke_size,
					  std::string stroke_color, std::string fill_color);

bool Write_SVG_rect(std::ostream& os, double UL_corner_x, double UL_corner_y, double width, double height,
					double stroke_size, std::string stroke_color, std::string fill_color);

bool Write_SVG_line(std::ostream& os, double start_x, double start_y, double end_x, double end_y, double thickness,
					std::string stroke_color);


} // namespace physicell
