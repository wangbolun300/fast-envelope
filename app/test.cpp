#include<fastenvelope/csv_reader.h>
#include<fastenvelope/FastEnvelope.h>
#include <igl/Timer.h>
#include <fastenvelope/getRSS.hpp>
#include <fstream>
#include <fastenvelope/MeshIO.hpp>

using namespace fastEnvelope;
//void pure_our_method(string queryfile, string model, string resultfile, Scalar envsize, bool csv_model) {
//	std::cout << "running our method" << std::endl;
//	vector<int> outenvelope;
//	std::vector<Vector3> env_vertices, v;
//	std::vector<Vector3i> env_faces, f;
//	GEO::Mesh envmesh, mesh;
//
//	std::vector<std::array<Vector3, 3>> triangles;
//	if (csv_model) {
//		triangles = read_CSV_triangle(queryfile, outenvelope);
//	}
//	else {
//		bool ok1 = MeshIO::load_mesh(queryfile, v, f, mesh);
//		if (!ok1) {
//			std::cout << ("Unable to load query mesh") << std::endl;
//			return;
//		}
//		triangles.resize(f.size());
//		for (int i = 0; i < f.size(); i++) {
//			for (int j = 0; j < 3; j++) {
//				triangles[i][j][0] = v[f[i][j]][0];
//				triangles[i][j][1] = v[f[i][j]][1];
//				triangles[i][j][2] = v[f[i][j]][2];
//			}
//		}
//	}
//
//	bool ok = MeshIO::load_mesh(model, env_vertices, env_faces, envmesh);
//	if (!ok) {
//		std::cout << ("Unable to load mesh") << std::endl;
//		return;
//	}
//
//	Vector3 min, max;
//	Scalar eps = envsize;
//
//	Scalar dd;
//	algorithms::get_bb_corners(env_vertices, min, max);
//	dd = ((max - min).norm()) *eps;
//	igl::Timer timer;
//	timer.start();
//	const FastEnvelope fast_envelope(env_vertices, env_faces, dd);
//
//	timer.stop();
//	const auto init_time = timer.getElapsedTimeInSec();
//	std::cout << "ours initialization time " << timer.getElapsedTimeInSec() << std::endl;
//	//int fn = triangles.size() > 100000 ? 100000 : triangles.size();
//	int fn = triangles.size() > 100000 ? 100000 : triangles.size();
//	std::cout << "total query size, " << fn << std::endl;
//	std::vector<bool> results;
//	results.resize(fn);
//	timer.start();
//	int inbr = 0;
//	for (int i = 0; i < fn; i++) {
//
//		results[i] = fast_envelope.is_outside(triangles[i]);
//		if (results[i] == 0) inbr++;
//	}
//	timer.stop();
//	cout << "ours method time " << timer.getElapsedTimeInSec() << endl;
//	cout << "memory use, " << getPeakRSS() << std::endl;
//	cout << "inside percentage, " << float(inbr) / float(fn) << std::endl;
//
//	std::ofstream fout;
//	fout.open(resultfile + ".json");
//	fout << "{\n";
//
//	fout << "\"method\": " << "\"ours\"" << ",\n";
//	fout << "\"init_time\": " << init_time << ",\n";
//	fout << "\"query_time\": " << timer.getElapsedTimeInSec() << ",\n";
//	fout << "\"memory\": " << getPeakRSS() << ",\n";
//	fout << "\"inside\": " << double(inbr) / double(fn) << ",\n";
//	fout << "\"queries\": " << fn << ",\n";
//	fout << "\"vertices\": " << env_vertices.size() << ",\n";
//	fout << "\"facets\": " << env_faces.size() << "\n";
//	fout << "}";
//	fout.close();
//	fast_envelope.printnumber();
//
//
//
//
//	fout.open(resultfile);
//	fout << "results" << endl;
//	for (int i = 0; i < fn; i++) {
//
//		fout << results[i] << endl;
//
//	}
//	fout.close();
//	std::cout << model << " done! " << std::endl;
//}

int main() {
	return 0;
}