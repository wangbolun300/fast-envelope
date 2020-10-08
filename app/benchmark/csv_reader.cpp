#include <fastenvelope/csv_reader.h>
#include<iostream>
#include<fstream>
#include<fastenvelope/mesh_AABB.h>
using namespace std;
namespace fastEnvelope{
	std::vector<std::array<Vector3, 3>> read_CSV_triangle(const string inputFileName, vector<int>& inenvelope) {


		std::vector<std::array<Vector3, 3>> triangle;

		ifstream infile;
		infile.open(inputFileName);
		if (!infile.is_open())
		{
			cout << "Path Wrong!!!!" << endl;
			return triangle;
		}

		int l = 0;
		while (infile) // there is input overload classfile
		{
			l++;
			string s;
			if (!getline(infile, s)) break;
			if (s[0] != '#') {
				istringstream ss(s);
				array<double, 10> record;
				int c = 0;
				while (ss) {
					string line;
					if (!getline(ss, line, ','))
						break;
					try {
						record[c] = stod(line);
						c++;
					}
					catch (const std::invalid_argument e) {
						cout << "NaN found in file " << inputFileName << " line " << l
							<< endl;
						e.what();
					}
				}

				triangle.push_back({ {Vector3(record[0],record[1],record[2]),Vector3(record[3],record[4],record[5]),
					Vector3(record[6],record[7],record[8])} });
				inenvelope.push_back(record[9]);
			}
		}
		if (!infile.eof()) {
			cerr << "Could not read file " << inputFileName << "\n";
		}
		cout << "triangle size " << triangle.size() << endl;
		return triangle;
	}

}