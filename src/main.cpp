#include <fastenvelope/FastEnvelope.h>
#include<iostream>
#include<array>
#include <fastenvelope/MeshIO.hpp>
#include <fstream>
#include <istream>
#include<fastenvelope/Predicates.hpp>
#include<fastenvelope/Types.hpp>
#include <fastenvelope/AABBWrapper.h>
#include <igl/Timer.h>
#include <ctime>
#include <cstdlib>
#include<fastenvelope/EnvelopeTest.h>
#include <unordered_map>

#include <filesystem>
#include <stdio.h>



using namespace floatTetWild;
using namespace fastEnvelope;
using namespace std;



void tri_tri_cutting_try() {

	Vector3 p1, p2, p3, q1, q2, q3;

	p1 = { 1,0,0 };

	p2 = { 0,1,0 };

	p3 = { 0,0,0 };


	q1 = { 0,0,0 };

	q2 = { 0,0.5,-1 };

	q3 = { 0,1,1 };

	int a = 0; //FIXME maybe
	//FastEnvelope::tri_cut_tri_simple(p1, p2, p3, q1, q3, q2);
	std::cout << "a " << a << std::endl;

}


void unordered_map_try() {
	std::vector<int> a, b;
	a.push_back(0);
	a.push_back(2);
	b.push_back(7);
	std::unordered_map<int, std::vector<int>> letter;
	letter[1] = a;
	letter[3] = b;
	letter[4].push_back(1);
	auto search = letter.find(4);
	if (search == letter.end()) {
		std::cout << "no find" << std::endl;
	}
	std::cout << "Found or not " << search->first << " size " << search->second.size() << '\n';
}
void get_bb_corners(const Parameters &params, const std::vector<Vector3> &vertices, Vector3 &min, Vector3 &max) {
	min = vertices.front();
	max = vertices.front();

	for (size_t j = 0; j < vertices.size(); j++) {
		for (int i = 0; i < 3; i++) {
			min(i) = std::min(min(i), vertices[j](i));
			max(i) = std::max(max(i), vertices[j](i));
		}
	}

	const Scalar dis = (max - min).minCoeff() * params.box_scale;

	for (int j = 0; j < 3; j++) {
		min[j] -= dis;
		max[j] += dis;
	}

	//cout << "min = " << min[0] << " " << min[1] << " " << min[2] << endl;
	//cout << "max = " << max[0] << " " << max[1] << " " << max[2] << endl;
	//            pausee();
}
void get_triangle_corners(const std::array<Vector3,3> &triangle, Vector3 &min, Vector3 &max) {
	min[0] = std::min(std::min(triangle[0][0], triangle[1][0]), triangle[2][0]);
	min[1] = std::min(std::min(triangle[0][1], triangle[1][1]), triangle[2][1]);
	min[2] = std::min(std::min(triangle[0][2], triangle[1][2]), triangle[2][2]);
	max[0] = std::max(std::max(triangle[0][0], triangle[1][0]), triangle[2][0]);
	max[1] = std::max(std::max(triangle[0][1], triangle[1][1]), triangle[2][1]);
	max[2] = std::max(std::max(triangle[0][2], triangle[1][2]), triangle[2][2]);

}
void CornerList(const Parameters& params, const std::vector<std::array<Vector3, 12>>& prism,
	std::vector<std::array<Vector3, 2>>& list) {
	std::vector<Vector3> ver12(12);
	Vector3 min, max;
	list.resize(prism.size());//to be safer
	for (int i = 0; i < prism.size(); i++) {
		for (int j = 0; j < 12; j++) {
			ver12[j] = prism[i][j];
		}
		get_bb_corners(params, ver12, min, max);
		list[i] = { {min,max } };
	}
}

void BondingBoxIntersectionFinding(Parameters& params, const std::array<Vector3, 3>& triangle,
	const std::vector<std::array<Vector3, 2>>& list, std::vector<int>& inumber) {
	inumber.clear();
	std::vector<Vector3> ver(3);
	inumber.reserve(list.size());
	Vector3 tmin, tmax;
	int rcd = 0;
	ver[0] = triangle[0];
	ver[1] = triangle[1];
	ver[2] = triangle[2];
	get_bb_corners(params, ver, tmin, tmax);
	for (int i = 0; i < list.size(); i++) {
		rcd = 0;
		for (int j = 0; j < 3; j++) {
			if (tmax[j] <= list[i][0][j] || tmin[j] >= list[i][1][j]) {
				rcd = 1;
				break;
			}

		}
		if (rcd == 0) {
			inumber.push_back(i);
		}
	}
}

//void BoxFindCells(const Vector3& min,const Vector3& max,
//	const Vector3& cellmin, const Vector3& cellmax, const int& sub, std::vector<int>& intercell) {
//	
//	Vector3 delta= (cellmax - cellmin) / sub;
//	//intercell.reserve(int((max - min)[0] / delta[0])*int((max - min)[1] / delta[1])*int((max - min)[2] / delta[2]));
//	intercell.clear();
//	int location[2][3];
//	for (int i = 0; i < 3; i++) {
//		location[0][i] = (min[i] - cellmin[i]) / delta[i];
//	}
//	for (int i = 0; i < 3; i++) {
//		location[1][i] = (max[i] - cellmin[i]) / delta[i];
//	}
//	for (int i = location[0][0]; i <= location[1][0]; i++) {
//		for (int j = location[0][1]; j <= location[1][1]; j++) {
//			for (int k = location[0][2]; k <= location[1][2]; k++) {
//				intercell.emplace_back(i*sub*sub + j * sub + k);
//			}
//		}
//	}
//
//	
//}
void BoxFindCellsV1(const Vector3& min, const Vector3& max,
	const Vector3& cellmin, const Vector3& cellmax, const int& subx,const int&suby, const int subz, std::vector<int>& intercell) {

	Vector3 delta;
	delta[0]= (cellmax - cellmin)[0] / subx;
	delta[1] = (cellmax - cellmin)[1] / suby;
	delta[2] = (cellmax - cellmin)[2] / subz;
	//intercell.reserve(int((max - min)[0] / delta[0])*int((max - min)[1] / delta[1])*int((max - min)[2] / delta[2]));
	intercell.clear();
	int location[2][3];
	for (int i = 0; i < 3; i++) {
		location[0][i] = (min[i] - cellmin[i]) / delta[i];
	}
	for (int i = 0; i < 3; i++) {
		location[1][i] = (max[i] - cellmin[i]) / delta[i];
	}
	for (int i = location[0][0]; i <= location[1][0]; i++) {
		for (int j = location[0][1]; j <= location[1][1]; j++) {
			for (int k = location[0][2]; k <= location[1][2]; k++) {
				intercell.emplace_back(k*subx*suby+j*subx+i);
			}
		}
	}


}


void testOrientation() {
	std::array<Vector3, 3> tri = {{ Vector3(0,0,0),Vector3(1,0,0),Vector3(0,1,0) }};// right hand law
	Vector3 p = { 0.4,0.4,-1 };
	int ori = Predicates::orient_3d(tri[0], tri[1], tri[2], p);
	std::cout << "orientation test : " << ori << std::endl;
}

bool is_out_function(const std::array<Vector3, 3>& triangle, const Scalar& dd, AABBWrapper& sf_tree) {
	std::vector<GEO::vec3> ps;
	sample_triangle(triangle, ps, dd);//dd is used for sapmling
	return sf_tree.is_out_sf_envelope(ps, pow(dd*(1 - (1 / sqrt(3))), 2));
	//int cnt = 0;
	//Scalar sq_dist;
	//double sq_distg;
	//GEO::vec3 nearest_point;
	//const unsigned int ps_size = ps.size();
	//for (unsigned int i = ps_size / 2;; i = (i + 1) % ps_size) {//check from the middle
	////for (unsigned int i = 0; i < ps.size(); i++) {
	//	GEO::vec3 &current_point = ps[i];
	//	sf_tree.nearest_facet(current_point, nearest_point, sq_distg);
	//	sq_dist = sq_distg;
	//	if (sq_dist > pow(dd*(1 - (1 / sqrt(3))), 2))// eps_here=eps-dd/sprt(3). in this case eps=dd
	//		return true;
	//	cnt++;

	//	if (cnt >= ps_size)
	//		break;


	//}
	//return false;
}






void add_hashing() {
	igl::Timer timer1, timer2, timer3, timer4, timer5;
	Scalar time1 = 0, time2 = 0, time3 = 0, time4 = 0, time5 = 0;
	const std::string input_surface_path1 = "D:\\vs\\float project\\data\\Roal10.stl";
	std::vector<Vector3> env_vertices;
	std::vector<Vector3i> env_faces;
	GEO::Mesh envmesh;
	
	///////////////////////////////////////////////////////
	bool ok1 = MeshIO::load_mesh(input_surface_path1, env_vertices, env_faces, envmesh);
	if (!ok1) {
		std::cout << ("Unable to load mesh") << std::endl;
		return;
	}
	std::cout << "envface size  " << env_faces.size() << "\nenv ver size " << env_vertices.size() << std::endl;
	/////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////
	//for aabb method
	Vector3 min, max;
	Parameters params;
	Scalar dd;
	get_bb_corners(params, env_vertices, min, max);
	///////////////////////////////////////////////////////

	Scalar shrink = 10;
	Scalar eps = 1e-3; 
	const int spac = 100;// space subdivision parameter
	const int fn = 80000;//test face number
	//////////////////////////////////////////////////////////////
	eps = eps / shrink; 
	eps=eps*sqrt(3)*(1 - (1 / sqrt(3)));//TODO to make bbd similar size to aabb method
	
	const FastEnvelope fast_envelope(env_vertices, env_faces, eps, spac);

	dd = ((max - min).norm()) / 1000 / shrink;
	std::vector<GEO::vec3> ps;
	int cnt = 0;
	double sq_distg;
	GEO::vec3 nearest_point;
	unsigned int ps_size;
	AABBWrapper sf_tree(envmesh);


	//////////////////////////generation of noised points
	/*
		std::ofstream fout;
	fout.open("D:\\vs\\float project\\data\\GenerateVer.obj");
	for (int i = 0; i < env_faces.size()*3; i++) {
		unsigned int time(0);
			fout <<"v " <<env_vertices[i][0] + dd * (rand() % 50-25)/shrink << " " << env_vertices[i][1]+ dd*(rand() % 50-25)/shrink << " " << env_vertices[i][2]+ dd*(rand() % 50-25)/shrink << endl;

	}
	for (int i = 0; i < env_faces.size(); i++) {

		fout << "f " << env_faces[i][0]+1 << " " << env_faces[i][1]+1 << " " << env_faces[i][2]+1 << endl;

	}
	fout.close();
	*/



	//////////////////////////generation of noised points finished


	const std::string input_surface_path2 = "D:\\vs\\float project\\data\\GenerateVer.obj";
	std::vector<Vector3> testvertices;
	std::vector<Vector3i> testfaces;
	GEO::Mesh testmesh;
	bool ok2 = MeshIO::load_mesh(input_surface_path2, testvertices, testfaces, testmesh);
	if (!ok2) {
		std::cout<<("Unable to load mesh")<<std::endl;
			return;
	}

	/////////////////////////////////

	bool pos1[fn], pos2[fn];
	std::vector<std::array<Vector3, 3>> triangle(fn);

	for (int i = 0; i < fn; i++) {
		triangle[i] = { testvertices[testfaces[i][0]],testvertices[testfaces[i][1]],    testvertices[testfaces[i][2]] };
	}
	/////////////////////////////////////////
	timer1.start();
	for (int i = 0; i < fn; i++) {
		
		pos1[i] = is_out_function(triangle[i], dd, sf_tree);	;
	}
	time1 = time1 + timer1.getElapsedTimeInSec();
	std::cout << "TEST ONE FINISHED  " << std::endl;
	//////////////////////////////
	



	timer2.start();

	for (int i = 0; i < fn; i++) {



		timer4.start();//function time

		pos2[i] = fast_envelope.is_outside(triangle[i]);
		time5 = time5 + timer4.getElapsedTimeInSec();//function time


	}
	time2 = time2 + timer2.getElapsedTimeInSec();

	std::cout << "TEST TWO FINISHED  " << std::endl;
	/////////////////////////////////
	/////////////////////////////////

	int rcd = 0,eq0=0,eq02=0;
	for (int i = 0; i < fn; i++) {
		if (pos1[i] - pos2[i] != 0) {
			//if (pos1[i]== 0) {
			rcd = rcd + 1;
			std::cout << "envelope test different! different face NO. " << i << " the difference: " << pos1[i] - pos2[i] << std::endl;
			//std::cout << "envelope test same! same face NO. " << i << "the in and out : " <<pos1[i] <<","<<pos2[i] << std::endl;
		}
		if (pos1[i] == 0) {
			eq0 = eq0 + 1;
		}
		if (pos2[i] == 0) {
			eq02 = eq02 + 1;
		}
	}
	std::cout << "aabb inside triangle number:  " << eq0 << std::endl;
	std::cout << "our  inside triangle number:  " << eq02 << std::endl;
	std::cout << "how many different cases:  " << rcd << std::endl;
	std::cout << "time1 and time2:  " << time1 << "," << time2 << std::endl;
	std::cout << "dd:  " << dd << std::endl;
	std::cout << "shrink size:  " << shrink << std::endl;
	std::cout << "all the prism size:  " << fast_envelope.prism_size() << std::endl;
	std::cout << "\ngetbbcorners time:  " << time3 << std::endl;
	std::cout << "intersection element finding time:  " << time4 << std::endl;
	std::cout << "function time:  " << time5 << std::endl;
	//FastEnvelope::timerecord();
	
	//13962,shrink 10

	///////////////
}
void calculation() {
	Eigen::Matrix<Scalar,4, 4> biga,bigb;
	Eigen::Matrix<Scalar, 3, 3> a;
	Eigen::Matrix<Scalar, 3, 4> b;
	Eigen::Matrix<Scalar, 3, 1> c0;
	Eigen::Matrix<Scalar, 1, 3> r0;
	Eigen::Matrix<Scalar, 4, 6> bigc;
	a << 1, 2, 3,
		4, 5, 6,
		7, 8, 9;
	c0 << 0, 
		0, 
		0;
	r0 << 0, 0, 0;
	biga << a, c0,
		r0, 1;
	b << 1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12;
	bigb << b,
		1, 1, 1, 1;
	//std::cout << biga*bigb<<"\n" << std::endl;
	//std::cout << a*b<< "\n" << std::endl;
	
	a.setZero();
	std::array<Vector3, 3> triangle0;
	triangle0[0] = { {Vector3(0,0,1)} };
	triangle0[1] = { {Vector3(0,1,1)} };
	triangle0[2] = { {Vector3(7,0,1)} };
	a.transpose() << (triangle0[0] - triangle0[1]).cross(triangle0[0] - triangle0[2]),
		(triangle0[0] - triangle0[1]).cross(triangle0[0] - triangle0[2]),
		(triangle0[0] - triangle0[1]).cross(triangle0[0] - triangle0[2]);
	bigc << a, a,
		1, 1, 1, 1, 1, 1;
	std::cout << a  << "\n" << std::endl;
	std::cout << bigc << "\n" << std::endl;
}

void test_ttt() {
	std::array<Vector3, 3> triangle,facet1,facet2,facet3;
	Eigen::Matrix<Scalar, 3, 3> A, AT, ATA;
	Eigen::Matrix<Scalar, 3, 1> B;
	Eigen::Matrix<Scalar, 4, 4> C;
	triangle = { { Vector3(0,0,0),Vector3(1,0,0),Vector3(0,1,0)} };
	facet1 = { { Vector3(0,0,0),Vector3(0,0,1),Vector3(1,0,0)} };
	facet2 = { { Vector3(0,0,0),Vector3(0,0,1),Vector3(0,1,0)} };
	facet3 = { { Vector3(0,0,0),Vector3(1,0,0),Vector3(0,1,0)} };
	AT << (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]),
		(facet1[0] - facet1[1]).cross(facet1[0] - facet1[2]),
		(facet2[0] - facet2[1]).cross(facet2[0] - facet2[2]);
	B << triangle[0], facet1[0], facet2[0];
	A = AT.transpose();
	B = A * B;
	ATA = AT * A;
	std::cout << A << std::endl;
	// if (FastEnvelope::determinant(A) == 0) { //mmnnmnmmnnnnnnnnnnnnnnnnnnnnnnmnmmmnmnmmmnmnmmmmmmnmn
	// 	std::cout << "A singular" << std::endl;
	// }
	int ori = -1; // FastEnvelope::orient_3triangles(A, AT, ATA, B, facet3); //FIXME maybe
	std::cout << ori << std::endl;
}
/*
void test_diff() {
	Vector3 fl1;
	std::vector<std::array<Vector3, 2>>seg;
	std::vector<Vector3> segun, triangleun, facetun;
	std::vector<Vector2i> segin;
	std::ifstream infile;
	Vector2i a;
	infile.open("D:\\vs\\float project\\data\\output_0522\\segpoints.txt");

	if (!infile.is_open())

		cout << "Open file failure" << endl;
	int k = 0;
	while (!infile.eof())

	{

		infile >> fl1[0] >> fl1[1] >> fl1[2];

		segun.push_back(fl1);
	}
	infile.close();
	std::cout << segun.size() << std::endl;
	std::cout << segun[2001] << std::endl;
	//std::cout << segun.size() << std::endl;

	/*infile.open("D:\\vs\\float project\\data\\output_0522\\segindex.txt");
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	k = 0;
	while (!infile.eof())

	{
		std::cout << k << std::endl;
		infile >> a[0] >> a[1];

		segin.push_back(a);
		k = k + 1;
	}
	infile.close();
	std::cout <<"segin size "<< segin.size() << std::endl;*/

	/*for (int i = 0; i < 1002; i++) {
		std::cout << i << " st\n " << segin[i] << std::endl;
	}

	infile.open("D:\\vs\\float project\\data\\output_0522\\triangle_points.txt");

	if (!infile.is_open())

		cout << "Open file failure" << endl;
	while (!infile.eof())

	{

		infile >> fl1[0] >> fl1[1] >> fl1[2];

		triangleun.push_back(fl1);
	}
	infile.close();
	cout << triangleun.size() << endl;



	infile.open("D:\\vs\\float project\\data\\output_0522\\facet_points.txt");

	if (!infile.is_open())

		cout << "Open file failure" << endl;
	while (!infile.eof())

	{

		infile >> fl1[0] >> fl1[1] >> fl1[2];

		facetun.push_back(fl1);
	}
	infile.close();
	cout << facetun.size() << endl;

	int num = segun.size() / 2;
	cout << num << endl;
	/////////////////////////////////////////////////////////////////////////////

	std::vector<std::array<Vector3, 2>> seglist;
	std::vector<std::array<Vector3, 3>> triangle, facet;
	for (int i = 0; i < num; i++) {
		seglist.push_back({ {segun[2 * i],segun[2 * i + 1]} });
		triangle.push_back({ {triangleun[3 * i],triangleun[3 * i + 1],triangleun[3 * i + 2]} });
		facet.push_back({ {facetun[3 * i],facetun[3 * i + 1],facetun[3 * i + 2]} });
	}
	cout << seglist.size() << " " << triangle.size() << " " << facet.size() << endl;
	//////////////////////////////////////////////////////////////////////////////////
	Scalar a11, a12, a13, a21, a22, a23, a31, a32, a33, px_rx, py_ry, pz_rz, d, n;
	for (int i = 0; i < num; i++) {
		bool inter = FastEnvelope::is_seg_facet_intersection(seglist[i][0][0], seglist[i][0][1], seglist[i][0][2],
			seglist[i][1][0], seglist[i][1][1], seglist[i][1][2],
			triangle[i][0][0], triangle[i][0][1], triangle[i][0][2],
			triangle[i][1][0], triangle[i][1][1], triangle[i][1][2],
			triangle[i][2][0], triangle[i][2][1], triangle[i][2][2],
			a11, a12, a13, a21, a22, a23, a31, a32, a33, px_rx, py_ry, pz_rz, d, n);
		if (inter == 0) {
			std::cout << "wrong in intersection" << endl;
		}
		int ori = FastEnvelope::orient3D_LPI(
			seglist[i][0][0], seglist[i][0][1], seglist[i][0][2],
			facet[i][0][0], facet[i][0][1], facet[i][0][2],
			facet[i][1][0], facet[i][1][1], facet[i][1][2],
			facet[i][2][0], facet[i][2][1], facet[i][2][2],
			a11, a12, a13, a21, a22, a23, a31, a32, a33, px_rx, py_ry, pz_rz, d, n);
		Vector3 point = seglist[i][0] + (n / d)*(seglist[i][0] - seglist[i][1]);
		int ori1 = -1 * Predicates::orient_3d(facet[i][0], facet[i][1], facet[i][2], point);
		if (ori == ori1) {
			std::cout << "out put ori wrong" << std::endl;
		}
		//std::cout << ori<<" "<<ori1 << endl;

	}





}
*/



std::vector<std::array<Vector3, 3>> read_CSV_triangle(const string inputFileName,vector<int>& inenvelope) {

	
	std::vector<std::array<Vector3,3>> triangle;
	
	ifstream infile;
	infile.open(inputFileName);
	if (!infile.is_open())
	{
		cout << "Path Wrong!!!!" << endl;
		exit(EXIT_FAILURE);
	}

	int l = 0;
	while (infile) // there is input overload classfile
	{
		l++;
		string s;
		if (!getline(infile, s)) break;
		if (s[0] != '#') {
			istringstream ss(s);
			array<double,10> record;
			int c = 0;
			while (ss) {
				string line;
				if (!getline(ss, line, ','))
					break;
				try {
					record[c]=stof(line);
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

void test_in_wild() {
	string inputFileName = "D:\\vs\\fast_envelope_csv\\thingi10k_debug\\100029\\100029.stl_env.csv";
	string input_surface_path1 = "D:\\vs\\fast_envelope_csv\\thingi10k_debug\\100029\\elevator_and_stabiliser_-_V4.stl";
	vector<int> outenvelope;
	std::vector<std::array<Vector3, 3>> triangles = read_CSV_triangle(inputFileName, outenvelope);
	
	std::vector<Vector3> env_vertices;
	std::vector<Vector3i> env_faces;
	GEO::Mesh envmesh;

	///////////////////////////////////////////////////////
	bool ok1 = MeshIO::load_mesh(input_surface_path1, env_vertices, env_faces, envmesh);
	if (!ok1) {
		std::cout << ("Unable to load mesh") << std::endl;
		return;
	}
	std::cout << "envface size  " << env_faces.size() << "\nenv ver size " << env_vertices.size() << std::endl;



	Scalar shrink = 10;
	Scalar eps = 1e-3;
	const int spac = 10;// space subdivision parameter
	const int fn = triangles.size();//test face number
	//////////////////////////////////////////////////////////////
	eps = eps / shrink;
	eps = eps * sqrt(3)*(1 - (1 / sqrt(3)));//TODO to make bbd similar size to aabb method

	const FastEnvelope fast_envelope(env_vertices, env_faces, eps, spac);
	vector<bool> pos1, pos2;
	pos1.resize(fn);
	pos2.resize(fn);
	for (int i = 0; i < fn; i++) {

		pos1[i] = outenvelope[i];
		pos2[i] = fast_envelope.is_outside(triangles[i]);

	}
	




	int rcd = 0, eq0 = 0, eq02 = 0,rmk=0;
	for (int i = 0; i < fn; i++) {
		if (pos1[i] - pos2[i] != 0) {
			//if (pos1[i]== 0) {
			rcd = rcd + 1;
			//std::cout << "envelope test different! different face NO. " << i << " the difference: " << pos1[i] - pos2[i] << std::endl;
			//std::cout << "envelope test same! same face NO. " << i << "the in and out : " <<pos1[i] <<","<<pos2[i] << std::endl;
		}
		if (pos1[i] == 0) {
			eq0 = eq0 + 1;
		}
		if (pos2[i] == 0) {
			eq02 = eq02 + 1;
		}
		if (pos1[i] == 0 && pos2[i] == 1) {
			rmk++;
		}
	}

	std::cout << "how many different cases:  " << rcd << std::endl;
	std::cout << "aabb inside triangle number:  " << eq0 << std::endl;
	std::cout << "our  inside triangle number:  " << eq02 << std::endl;
	std::cout << "0-1 cases number " << rmk << std::endl;



	 
	
	std::vector<bool> pos3;
	pos3.resize(fn);
	std::array<Vector3, 3> tri;
	Scalar l1;
	std::vector<std::array<Vector3, 3>> t;
	std::vector<int> trindex1, trindex2;

	for (int i = 0; i < fn; i++) {
		tri = triangles[i];
		l1 = (tri[0] - tri[1]).norm() / 30;
		pos3[i]=fast_envelope.sample_triangle_outside(tri, l1);
	}

	int count = 0, count1 = 0, count2 = 0, count3 = 0;
	for (int i = 0; i < fn; i++) {
		if (pos2[i] - pos3[i] != 0) {
			count++;
		}
		if (pos2[i] - pos3[i] == 1) {
			count1++;
		}
		if (pos2[i] - pos3[i] == 1) {
			trindex1.push_back(i);
		}
	}
	std::cout << "how many different cases in comparison:  " << count << std::endl;
	std::cout << "1-0 cases in comparison:  " << count1 << std::endl;
	std::cout << "0-1 cases in comparison:  " << count-count1 << std::endl;
	std::cout << "0-1 case size:  " << trindex1.size() << std::endl;

	for (int i = 0; i < fn; i++) {
		if (pos1[i] - pos3[i] != 0) {
			count2++;
		}
		if (pos1[i] - pos3[i] == 1) {
			count3++;
		}
	}
	std::cout << "how many different cases in comparison 1:  " << count2 << std::endl;
	std::cout << "1-0 cases in comparison 1:  " << count3 << std::endl;
	std::cout << "0-1 cases in comparison 1:  " << count2 - count3 << std::endl;


	////////////////////////////////////////////////////////////////////////////////////////
	std::ofstream fout;
	fout.open("D:\\vs\\fast_envelope_csv\\thingi10k_debug\\100029\\visualtriangle.txt");
	int idx = 2;
	std::cout << "output NO. " <<trindex1[idx] << endl;
	for (int i = 0; i < 3; i++) {

		fout << std::setprecision(17) << triangles[trindex1[idx]][i][0] << " " << triangles[trindex1[idx]][i][1] << " " << triangles[trindex1[idx]][i][2] << endl;

	}
	fout.close();

	fast_envelope.print_prisms(triangles[trindex1[idx]]);

	////////////////////////////////////////////////////////////////////////////////////////
}


void sample_triangle_test() {
	string inputFileName = "D:\\vs\\fast_envelope_csv\\thingi10k_debug\\100029\\100029.stl_env.csv";
	string input_surface_path1 = "D:\\vs\\fast_envelope_csv\\thingi10k_debug\\100029\\elevator_and_stabiliser_-_V4.stl";
	vector<int> outenvelope;
	std::vector<std::array<Vector3, 3>> triangles = read_CSV_triangle(inputFileName, outenvelope);
	std::array<Vector3, 3> tri = triangles[10000];
	std::vector<Vector3> ps;
	Scalar l1 = (tri[0] - tri[1]).norm()/50;
	FastEnvelope::triangle_sample(tri, ps, l1);
	std::cout << ps.size() << endl;


	std::ofstream fout;
	fout.open("D:\\vs\\fast_envelope_csv\\thingi10k_debug\\100029\\triangle.txt");
	for (int i = 0; i < 3; i++) {
		
		fout <<std::setprecision(17) << tri[i][0]<<" "<< tri[i][1]<<" "<< tri[i][2] <<endl;

	}
	fout.close();

	fout.open("D:\\vs\\fast_envelope_csv\\thingi10k_debug\\100029\\points.txt");
	for (int i = 0; i < ps.size(); i++) {

		fout << std::setprecision(17)<< ps[i][0] << " " << ps[i][1] << " " << ps[i][2] << endl;

	}
	fout.close();

}

//#include <fastenvelope/fprecision.h>

int main(int argc, char const *argv[])
{
	GEO::initialize();

	/*int n_digits = 400;

	float_precision num1(1, n_digits);
	float_precision num2(3, n_digits);

	float_precision res = num1/num2*float_precision(2) - float_precision(2,n_digits)/float_precision(3,n_digits);

	std::cout<<res<<std::endl;*/


	//assert(false);
	
	//EnvelopPrism();
	//PrismCutTriangle();
	//bool b= booltest();
	//std::cout << "bool test :" << b << std::endl;
	//TriangleInEnvelope();
	//RealModelEnvelopeTest();
	//testOrientation();
	//bunnytest();
	//change_model();
	//orient3D_LPITest();
	//change_model_with_prediction();
	//EnvelopeWithTree();
	//comparison();
	//unordered_map_try();
	//add_hashing();
	//tri_tri_cutting_try();
	//FastEnvelope::timerecord();
	//calculation();
	//test_ttt();
	//test_diff();
	test_in_wild();
	//inf();
	//sample_triangle_test();
	
	std::cout << "done!" << std::endl;
	

	return 0;
}
