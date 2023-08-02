#include "vector_utils.h"

using namespace biofvm;

namespace physicell {

std::vector<real_t> csv_to_vector(const std::string& buffer)
{
	std::vector<real_t> vect;
	unsigned int i = 0;
	while (i < buffer.size())
	{
		// churn through delimiters, whitespace, etc. to reach the next numeric term
		while (isdigit(buffer[i]) == false && buffer[i] != '.' && buffer[i] != '-' && buffer[i] != 'e'
			   && buffer[i] != 'E')
		{
			i++;
		}

		if (i < buffer.size()) // add this extra check in case of a final character, e.g., ']'
		{
			char* pEnd;
			vect.push_back(strtod(buffer.data() + i, &pEnd));
			i = pEnd - buffer.data();
		}
	}
	return vect;
}

void data_to_list(biofvm::real_t* data, std::size_t count, char*& buffer, char delim)
{
	// %.7e is approximately the same at matlab longe for single precision.
	// If you want better precision, use a binary data format like matlab, or (in the future) HDF

	int position = 0;
	for (unsigned int j = 0; j < count - 1; j++)
	{
		position += std::sprintf(buffer + position, "%.7e%c", (double)data[j], delim);
	}
	std::sprintf(buffer + position, "%.7e", data[count - 1]);
	return;
}

} // namespace physicell