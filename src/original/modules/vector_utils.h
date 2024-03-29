#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <BioFVM/types.h>

namespace physicell {

std::vector<biofvm::real_t> csv_to_vector(const std::string& value);

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
	for (unsigned int i = 0; i < v.size(); i++)
		os << v[i] << " ";

	return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const biofvm::point_t<T, 3>& v)
{
	for (unsigned int i = 0; i < v.size(); i++)
		os << v[i] << " ";

	return os;
}

void data_to_list(biofvm::real_t* data, std::size_t count, char*& buffer, char delim);

biofvm::real_t norm(const biofvm::point_t<biofvm::real_t, 3>& vec);

} // namespace physicell
