#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/MeshIO.hpp>


#include <fastenvelope/Types.hpp>
#include <fastenvelope/AABBWrapper.h>


#include <igl/Timer.h>


#include <geogram/basic/command_line.h>

#include <geogram/basic/command_line_args.h>

#ifdef USE_TBB
#include <tbb/tbb.h>
#endif

#include <unordered_map>

#include <stdlib.h>

#include <stdio.h>
#include <assert.h>
#include <fstream>
#include <istream>
#include <iostream>
#include <array>
#include <ctime>
#include <cstdlib>
#include<igl/writeOBJ.h>
#include<igl/readSTL.h>
#include<igl/writeOFF.h>
#include <fastenvelope/getRSS.hpp>

#ifdef ENVELOPE_WITH_GMP
#include <fastenvelope/Rational.hpp>
#endif
using namespace fastEnvelope;
using namespace std;



void get_bb_corners(const std::vector<Vector3> &vertices, Vector3 &min, Vector3 &max) {
	min = vertices.front();
	max = vertices.front();

	for (size_t j = 0; j < vertices.size(); j++) {
		for (int i = 0; i < 3; i++) {
			min(i) = std::min(min(i), vertices[j](i));
			max(i) = std::max(max(i), vertices[j](i));
		}
	}
}

bool sample_trianglex(const std::array<Vector3, 3>& vs, std::vector<GEO::vec3>& ps, Scalar sampling_dist, AABBWrapper& sf_tree) {
    
    GEO::vec3 nearest_point;
    double sq_dist = std::numeric_limits<double>::max();
    GEO::index_t prev_facet = GEO::NO_FACET;
    Scalar eps_2 = pow(sampling_dist*(1 - (1 / sqrt(3))), 2);


    Scalar sqrt3_2 = std::sqrt(3) / 2;

	std::array<Scalar, 3> ls;
	for (int i = 0; i < 3; i++) {
		ls[i] = (vs[i] - vs[(i + 1) % 3]).squaredNorm();
	}
	auto min_max = std::minmax_element(ls.begin(), ls.end());
	int min_i = min_max.first - ls.begin();
	int max_i = min_max.second - ls.begin();
	Scalar N = sqrt(ls[max_i]) / sampling_dist;
	if (N <= 1) {
		for (int i = 0; i < 3; i++) {
//            ps.push_back(GEO::vec3(vs[i][0], vs[i][1], vs[i][2]));
            if(sf_tree.is_out_sf_envelope(vs[i], eps_2, prev_facet, sq_dist, nearest_point)) {
                return true;
            }
        }
		return false;
	}
	if (N == int(N))
		N -= 1;

	GEO::vec3 v0(vs[max_i][0], vs[max_i][1], vs[max_i][2]);
	GEO::vec3 v1(vs[(max_i + 1) % 3][0], vs[(max_i + 1) % 3][1], vs[(max_i + 1) % 3][2]);
	GEO::vec3 v2(vs[(max_i + 2) % 3][0], vs[(max_i + 2) % 3][1], vs[(max_i + 2) % 3][2]);

	GEO::vec3 n_v0v1 = GEO::normalize(v1 - v0);
	for (int n = 0; n <= N; n++) {
//		ps.push_back(v0 + n_v0v1 * sampling_dist * n);
        if(sf_tree.is_out_sf_envelope(v0 + n_v0v1 * sampling_dist * n, eps_2, prev_facet, sq_dist, nearest_point)) {
            return true;
        }
	}
//	ps.push_back(v1);
    if(sf_tree.is_out_sf_envelope(v1, eps_2, prev_facet, sq_dist, nearest_point)) {
        return true;
    }

	Scalar h = GEO::distance(GEO::dot((v2 - v0), (v1 - v0)) * (v1 - v0) / ls[max_i] + v0, v2);
	int M = h / (sqrt3_2 * sampling_dist);
	if (M < 1) {
		ps.push_back(v2);
		return false;
//        return sf_tree.is_out_sf_envelope(v2, eps_2, prev_facet, sq_dist, nearest_point);
	}

	GEO::vec3 n_v0v2 = GEO::normalize(v2 - v0);
	GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
	Scalar tan_v0, tan_v1, sin_v0, sin_v1;
	sin_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) / (GEO::distance(v0, v2) * GEO::distance(v0, v1));
	tan_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) / GEO::dot((v2 - v0), (v1 - v0));
	tan_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) / GEO::dot((v2 - v1), (v0 - v1));
	sin_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) / (GEO::distance(v1, v2) * GEO::distance(v0, v1));

	for (int m = 1; m <= M; m++) {
		int n = sqrt3_2 / tan_v0 * m + 0.5;
		int n1 = sqrt3_2 / tan_v0 * m;
		if (m % 2 == 0 && n == n1) {
			n += 1;
		}
		GEO::vec3 v0_m = v0 + m * sqrt3_2 * sampling_dist / sin_v0 * n_v0v2;
		GEO::vec3 v1_m = v1 + m * sqrt3_2 * sampling_dist / sin_v1 * n_v1v2;
		if (GEO::distance(v0_m, v1_m) <= sampling_dist)
			break;

		Scalar delta_d = ((n + (m % 2) / 2.0) - m * sqrt3_2 / tan_v0) * sampling_dist;
		GEO::vec3 v = v0_m + delta_d * n_v0v1;
		int N1 = GEO::distance(v, v1_m) / sampling_dist;
		for (int i = 0; i <= N1; i++) {
//			ps.push_back(v + i * n_v0v1 * sampling_dist);
            if(sf_tree.is_out_sf_envelope(v + i * n_v0v1 * sampling_dist, eps_2, prev_facet, sq_dist, nearest_point))
                return true;
		}
	}
//	ps.push_back(v2);
    if(sf_tree.is_out_sf_envelope(v2, eps_2, prev_facet, sq_dist, nearest_point))
        return true;

	//sample edges
	N = sqrt(ls[(max_i + 1) % 3]) / sampling_dist;
	if (N > 1) {
		if (N == int(N))
			N -= 1;
		GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
		for (int n = 1; n <= N; n++) {
//			ps.push_back(v1 + n_v1v2 * sampling_dist * n);
            if(sf_tree.is_out_sf_envelope(v1 + n_v1v2 * sampling_dist * n, eps_2, prev_facet, sq_dist, nearest_point))
                return true;
		}
	}

	N = sqrt(ls[(max_i + 2) % 3]) / sampling_dist;
	if (N > 1) {
		if (N == int(N))
			N -= 1;
		GEO::vec3 n_v2v0 = GEO::normalize(v0 - v2);
		for (int n = 1; n <= N; n++) {
//			ps.push_back(v2 + n_v2v0 * sampling_dist * n);
            if(sf_tree.is_out_sf_envelope(v2 + n_v2v0 * sampling_dist * n, eps_2, prev_facet, sq_dist, nearest_point))
                return true;
		}
	}

	return false;
}


bool is_out_function(const std::array<Vector3, 3>& triangle, const Scalar& dd, AABBWrapper& sf_tree) {
	std::vector<GEO::vec3> ps;
	return sample_trianglex(triangle, ps, dd, sf_tree);//dd is used for sapmling
}

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



//can get the time, memory and all the result for queries of sampling method
#include <geogram/mesh/mesh_AABB.h>
void pure_sampling(string queryfile, string model,string resultfile, Scalar envsize, bool csv_model) {
	std::cout << "using sampling method" << std::endl;
	vector<int> outenvelope;
	std::vector<Vector3> env_vertices, v;
	std::vector<Vector3i> env_faces, f;
	GEO::Mesh envmesh, mesh;


	


	std::vector<std::array<Vector3, 3>> triangles;
	if (csv_model) {
		triangles = read_CSV_triangle(queryfile, outenvelope);
	}
	else {
		bool ok1 = MeshIO::load_mesh(queryfile, v, f, mesh);
		if (!ok1) {
			std::cout << ("Unable to load query mesh") << std::endl;
			return;
		}
		triangles.resize(f.size());
		for (int i = 0; i < f.size(); i++) {
			for (int j = 0; j < 3; j++) {
				triangles[i][j] = v[f[i][j]];
			}
		}
	}

	bool ok = MeshIO::load_mesh(model, env_vertices, env_faces, envmesh);
	if (!ok) {
		std::cout << ("Unable to load mesh") << std::endl;
		return;
	}

	Vector3 min, max;
	Scalar eps = envsize;
	
	Scalar dd;
	get_bb_corners(env_vertices, min, max);
	dd = ((max - min).norm()) *eps;
	igl::Timer timer;
	timer.start();
	AABBWrapper sf_tree(envmesh);



	const double init_time = timer.getElapsedTimeInSec();
	std::cout << "sampling initialization time " << init_time<< std::endl;
	int fn = triangles.size() > 100000 ? 100000 : triangles.size();
	std::cout << "total query size, " << fn << std::endl;
	std::vector<bool> results;
	results.resize(fn);
	timer.start();
	int inbr = 0;
	cout<<"dd = "<<dd<<endl;
	for (int i = 0; i < fn; i++) {

		results[i] = is_out_function(triangles[i], dd, sf_tree); ;
		if (results[i] == 0) inbr++;

	}
	cout << "sampling method time " << timer.getElapsedTimeInSec() << endl;
	cout << "memory use, " << getPeakRSS() << std::endl;
	cout <<"inside% = "<<(inbr/double(fn))<<endl;
	std::ofstream fout;

	
	fout.open(resultfile + ".json");
	fout << "{\n";

	fout << "\"method\": " << "\"sampling\"" << ",\n";
	fout << "\"init_time\": " << init_time << ",\n";
	fout << "\"query_time\": " << timer.getElapsedTimeInSec() << ",\n";
	fout << "\"memory\": " << getPeakRSS() << ",\n";
	fout << "\"inside\": " << double(inbr) / double(fn) << ",\n";
	fout << "\"queries\": " << fn << ",\n";
	fout << "\"vertices\": " << env_vertices.size() << ",\n";
	fout << "\"facets\": " << env_faces.size() << "\n";
	fout << "}";
	fout.close();


	fout.open(resultfile);
	fout << "results" << endl;
	for (int i = 0; i < fn; i++) {

		fout << results[i] << endl;

	}
	fout.close();
	std::cout << model << " done! " << std::endl;
}

void pure_our_method(string queryfile, string model, string resultfile, Scalar envsize, bool csv_model) {
	std::cout << "running our method" << std::endl;
	vector<int> outenvelope;
	std::vector<Vector3> env_vertices, v;
	std::vector<Vector3i> env_faces, f;
	GEO::Mesh envmesh, mesh;

	std::vector<std::array<Vector3, 3>> triangles;
	if (csv_model) {
		triangles = read_CSV_triangle(queryfile, outenvelope);
	}
	else {
		bool ok1 = MeshIO::load_mesh(queryfile, v, f, mesh);
		if (!ok1) {
			std::cout << ("Unable to load query mesh") << std::endl;
			return;
		}
		triangles.resize(f.size());
		for (int i = 0; i < f.size(); i++) {
			for (int j = 0; j < 3; j++) {
				triangles[i][j][0] = v[f[i][j]][0];
				triangles[i][j][1] = v[f[i][j]][1];
				triangles[i][j][2] = v[f[i][j]][2];
			}
		}
	}

	bool ok = MeshIO::load_mesh(model, env_vertices, env_faces, envmesh);
	if (!ok) {
		std::cout << ("Unable to load mesh") << std::endl;
		return;
	}

	Vector3 min, max;
	Scalar eps = envsize;
	
	Scalar dd;
	get_bb_corners(env_vertices, min, max);
	dd = ((max - min).norm()) *eps;
	igl::Timer timer;
	timer.start();
	const FastEnvelope fast_envelope(env_vertices, env_faces, dd);

	timer.stop();
	const auto init_time = timer.getElapsedTimeInSec();
	std::cout << "ours initialization time " << timer.getElapsedTimeInSec() << std::endl;
	//int fn = triangles.size() > 100000 ? 100000 : triangles.size();
	int fn = triangles.size() > 1000 ? 1000 : triangles.size();
	std::cout << "total query size, " << fn << std::endl;
	std::vector<bool> results;
	results.resize(fn);
	timer.start();
	int inbr = 0;
	for (int i = 0; i < fn; i++) {

		results[i] = fast_envelope.is_outside(triangles[i]);
		if (results[i] == 0) inbr++;
	}
	timer.stop();
	cout << "ours method time " << timer.getElapsedTimeInSec() << endl;
	cout << "memory use, " << getPeakRSS() << std::endl;
	cout << "inside percentage, " << float(inbr) / float(fn) << std::endl;

	std::ofstream fout;
	fout.open(resultfile + ".json");
	fout << "{\n";

	fout << "\"method\": " << "\"ours\"" << ",\n";
	fout << "\"init_time\": " << init_time << ",\n";
	fout << "\"query_time\": " << timer.getElapsedTimeInSec() << ",\n";
	fout << "\"memory\": " << getPeakRSS() << ",\n";
	fout << "\"inside\": " << double(inbr) / double(fn) << ",\n";
	fout << "\"queries\": " << fn << ",\n";
	fout << "\"vertices\": " << env_vertices.size() << ",\n";
	fout << "\"facets\": " << env_faces.size() << "\n";
	fout << "}";
	fout.close();
	fast_envelope.printnumber();




	fout.open(resultfile);
	fout << "results" << endl;
	for (int i = 0; i < fn; i++) {

		fout << results[i] << endl;

	}
	fout.close();
	std::cout << model << " done! " << std::endl;
}


void read_CSV_triangle_write_csv(const string inputFileName) {

	std::vector<std::array<Vector3, 3>> triangle;
	


	ifstream infile;
	infile.open(inputFileName);
	if (!infile.is_open())
	{
		cout << "Path Wrong!!!!" << endl;
		return;
	}
	std::ofstream fout;
	fout.open("D:\\vs\\fast_envelope_csv\\python\\differenteps\\short.csv");

	int l = 0;
	while (infile) // there is input overload classfile
	{
		l++;
		string s;
		if (!getline(infile, s)) break;
		if (s[0] != '#') {
			istringstream ss(s);
			array<double, 11> record;
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
			if (record[10] < 1) {
				fout << std::setprecision(17) << record[0]<<","<<record[1] << ","<<record[2] << ","
					<<record[3] << "," <<record[4] << "," <<record[5] << "," 
					<<record[6] << "," <<record[7] << "," <<record[8] << "," <<record[9]  << endl;
			}
			
		}
	}
	if (!infile.eof()) {
		cerr << "Could not read file " << inputFileName << "\n";
	}
	cout << "triangle size " << triangle.size() << endl;
	fout.close();
	infile.close();
}
void pure_our_method_detailed(string queryfile, string model, string resultfile, Scalar envsize, bool csv_model) {
	std::cout << "running our method with detailed time information in csv" << std::endl;
	vector<int> outenvelope;
	std::vector<Vector3> env_vertices, v;
	std::vector<Vector3i> env_faces, f;
	GEO::Mesh envmesh, mesh;

	///////////////////////////////////////////////////////



	std::vector<std::array<Vector3, 3>> triangles;
	if (csv_model) {
		triangles = read_CSV_triangle(queryfile, outenvelope);
	}
	else {
		bool ok1 = MeshIO::load_mesh(queryfile, v, f, mesh);
		if (!ok1) {
			std::cout << ("Unable to load query mesh") << std::endl;
			return;
		}
		triangles.resize(f.size());
		for (int i = 0; i < f.size(); i++) {
			for (int j = 0; j < 3; j++) {
				triangles[i][j][0] = v[f[i][j]][0];
				triangles[i][j][1] = v[f[i][j]][1];
				triangles[i][j][2] = v[f[i][j]][2];
			}
		}
	}

	bool ok = MeshIO::load_mesh(model, env_vertices, env_faces, envmesh);
	if (!ok) {
		std::cout << ("Unable to load mesh") << std::endl;
		return;
	}

	Vector3 min, max;
	Scalar eps = envsize;
	
	Scalar dd;
	get_bb_corners(env_vertices, min, max);
	dd = ((max - min).norm()) *eps;
	igl::Timer timer, timerdetail;
	timer.start();
	/////////////////////////////////////////////////////////////
	FastEnvelope fast_envelope;
	fast_envelope.init(env_vertices, env_faces, dd);
	/////////////////////////////////////////////////////////////
	//FastEnvelope fast_envelope;
	//std::vector<Scalar> ddlist;
	//ddlist.resize(env_faces.size());
	//for (int i = 0; i < env_faces.size(); i++) {
	//	ddlist[i] = dd;
	//}
	//std::cout<<"using adaptive"<<std::endl;
	//fast_envelope.init(env_vertices, env_faces, ddlist);
	/////////////////////////////////////////////////////////////
	timer.stop();
	const auto init_time = timer.getElapsedTimeInSec();
	std::cout << "ours initialization time " << timer.getElapsedTimeInSec() << std::endl;
	int fn = triangles.size() > 100000 ? 100000 : triangles.size();
	std::cout << "total query size, " << fn << std::endl;
	std::vector<bool> results;
	std::vector<double> timerlist;
	timerlist.resize(fn);
	results.resize(fn);
	
	timer.start();
	int inbr = 0;
	for (int i = 0; i < fn; i++) {
		timerdetail.start();
		results[i] = 0;
		results[i] = fast_envelope.is_outside(triangles[i]);
		timerlist[i] = timerdetail.getElapsedTimeInSec();
		if (results[i] == 0) inbr++;
		//cout << "result for this one " << i << " " << results[i] << std::endl;
	}
	timer.stop();
	cout << "ours method time " << timer.getElapsedTimeInSec() << endl;
	cout << "memory use, " << getPeakRSS() << std::endl;
	cout << "inside percentage, " << float(inbr) / float(fn) << std::endl;

	std::ofstream fout;
	fout.open(resultfile + "_detail.csv");
	fout << "triangle,time, out" << std::endl;
	for (int i = 0; i < fn; i++) {
		fout << i << "," << timerlist[i] << "," << results[i] << std::endl;
	}

	fout.close();
}


void pure_our_method_no_optimization(string queryfile, string model, string resultfile, Scalar envsize, bool csv_model) {
	std::cout << "running our method without optimazation" << std::endl;
	vector<int> outenvelope;
	std::vector<Vector3> env_vertices, v;
	std::vector<Vector3i> env_faces, f;
	GEO::Mesh envmesh, mesh;

	///////////////////////////////////////////////////////



	std::vector<std::array<Vector3, 3>> triangles;
	if (csv_model) {
		triangles = read_CSV_triangle(queryfile, outenvelope);
	}
	else {
		bool ok1 = MeshIO::load_mesh(queryfile, v, f, mesh);
		if (!ok1) {
			std::cout << ("Unable to load query mesh") << std::endl;
			return;
		}
		triangles.resize(f.size());
		for (int i = 0; i < f.size(); i++) {
			for (int j = 0; j < 3; j++) {
				triangles[i][j][0] = v[f[i][j]][0];
				triangles[i][j][1] = v[f[i][j]][1];
				triangles[i][j][2] = v[f[i][j]][2];
			}
		}
	}

	bool ok = MeshIO::load_mesh(model, env_vertices, env_faces, envmesh);
	if (!ok) {
		std::cout << ("Unable to load mesh") << std::endl;
		return;
	}

	Vector3 min, max;
	Scalar eps = envsize;
	
	Scalar dd;
	get_bb_corners(env_vertices, min, max);
	dd = ((max - min).norm()) *eps;
	igl::Timer timer;
	timer.start();
	const FastEnvelope fast_envelope(env_vertices, env_faces, dd);
	
	timer.stop();
	const auto init_time = timer.getElapsedTimeInSec();
	std::cout << "ours initialization time " << timer.getElapsedTimeInSec() << std::endl;
	int fn = triangles.size() > 100000 ? 100000 : triangles.size();
	std::cout << "total query size, " << fn << std::endl;
	std::vector<bool> results;
	results.resize(fn);
	timer.start();
	int inbr = 0;
	for (int i = 0; i < fn; i++) {

		results[i] = fast_envelope.is_outside_no_optimazation(triangles[i]);
		if (results[i] == 0) inbr++;
	}
	timer.stop();
	cout << "ours method time " << timer.getElapsedTimeInSec() << endl;
	cout << "memory use, " << getPeakRSS() << std::endl;
	cout << "inside percentage, " << float(inbr) / float(fn) << std::endl;

	std::ofstream fout;
	fout.open(resultfile + ".json");
	fout << "{\n";

	fout << "\"method\": " << "\"ours\"" << ",\n";
	fout << "\"init_time\": " << init_time << ",\n";
	fout << "\"query_time\": " << timer.getElapsedTimeInSec() << ",\n";
	fout << "\"memory\": " << getPeakRSS() << ",\n";
	fout << "\"inside\": " << double(inbr) / double(fn) << ",\n";
	fout << "\"queries\": " << fn << ",\n";
	fout << "\"vertices\": " << env_vertices.size() << ",\n";
	fout << "\"facets\": " << env_faces.size() << "\n";
	fout << "}";
	fout.close();





	fout.open(resultfile);
	fout << "results" << endl;
	for (int i = 0; i < fn; i++) {

		fout << results[i] << endl;

	}
	fout.close();
	std::cout << model << " done! " << std::endl;
}
//void appendix() {
//	Vector3 s(1, 1, 1);
//	Vector3 r(2, 4, 3);
//	Vector3 t(0, 3, 7);
//	Vector3 n = (r - s).cross(r - t);
//	std::cout << "normal " << n << std::endl;
//	Scalar p1 = 1, p2 = 1, p3;
//	p3 = n.dot(r)-n[0]*p1-n[1]
//	Vector3 bary = (v1 + v2 + v3) / 3;
//	std::cout <<"ori " <<Predicates::orient_3d(bary, v1, v2, v3) << std::endl;
//	std::cout << std::setprecision(16)<< "bary coor " << bary[0] << " " << bary[1] << " " << bary[2] << std::endl;
//
//}
//

int main(int argc, char const *argv[])
{
	//srand(42);

#ifndef WIN32
	setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif
	GEO::initialize(0);
	GEO::CmdLine::import_arg_group("standard");
	GEO::CmdLine::import_arg_group("pre");
	GEO::CmdLine::import_arg_group("algo");
#ifdef ENVELOPE_WITH_GMP
	std::cout << "using RATIONAL calculation in GMP" << std::endl;
#endif

	//string queryfile, string model, string resultfile, Scalar envelope size ratio epsilon, bool csv_model
	string keyword = argv[6];
	if (keyword == "sampling") {
		pure_sampling(argv[1], argv[2], argv[3], stod(argv[4]), stoi(argv[5]));
	}
	else if (keyword == "ours") {
		pure_our_method(argv[1], argv[2], argv[3], stod(argv[4]), stoi(argv[5]));
	}
	else if (keyword == "ours_without_optimazation") {
		pure_our_method_no_optimization(argv[1], argv[2], argv[3], stod(argv[4]), stoi(argv[5]));
	}
	else if (keyword == "ours_detailed_time") {
		pure_our_method_detailed(argv[1], argv[2], argv[3], stod(argv[4]), stoi(argv[5]));
	}
	else {
		std::cout << "wrong arguments" << std::endl;
	}
	
	return 0;
}
