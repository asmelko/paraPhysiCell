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

} // namespace physicell