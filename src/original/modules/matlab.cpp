#define _CRT_SECURE_NO_WARNINGS

#include "matlab.h"

#include <iostream>

namespace physicell {

FILE* write_matlab4_header(int nrows, int ncols, std::string filename, std::string variable_name)
{
	FILE* fp;
	fp = fopen(filename.c_str(), "wb");
	if (fp == NULL)
	{
		std::cout << "Error: could not open file " << filename << "!" << std::endl;
		return NULL;
	}

	typedef unsigned int UINT;
	UINT UINTs = sizeof(UINT);

	UINT temp;

	UINT type_numeric_format = 0; // little-endian assumed for now!
	UINT type_reserved = 0;
	UINT type_data_format = 0; // doubles for all entries
	UINT type_matrix_type = 0; // full matrix, not sparse

	temp = 1000 * type_numeric_format + 100 * type_reserved + 10 * type_data_format + type_matrix_type;
	fwrite((char*)&temp, UINTs, 1, fp);

	// UINT rows = (UINT) number_of_data_entries; // storing data as rows
	UINT rows = (UINT)nrows; // size_of_each_datum; // storing data as cols
	fwrite((char*)&rows, UINTs, 1, fp);

	// UINT cols = (UINT) size_of_each_datum; // storing data as rows
	UINT cols = (UINT)ncols; // number_of_data_entries; // storing data as cols
	fwrite((char*)&cols, UINTs, 1, fp);

	UINT imag = 0; // no complex matrices!
	fwrite((char*)&imag, UINTs, 1, fp);

	UINT name_length = variable_name.size(); // strlen( variable_name );
	fwrite((char*)&name_length, UINTs, 1, fp);

	// this is the end of the 20-byte header

	// write the name

	fwrite(variable_name.c_str(), name_length, 1, fp);

	return fp;
}

FILE* write_matlab_header(unsigned int rows, unsigned int cols, std::string filename, std::string variable_name)
{
	return write_matlab4_header(rows, cols, filename, variable_name);
}

} // namespace physicell
