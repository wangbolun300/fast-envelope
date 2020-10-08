#pragma once
#include <fastenvelope/Types.hpp>
#include <vector>
#include <array>
using namespace std;
namespace fastEnvelope {
	std::vector<std::array<Vector3, 3>> read_CSV_triangle(const string inputFileName, vector<int>& inenvelope);
}