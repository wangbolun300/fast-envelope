#include<fastenvelope/obb.h>
#include <Eigen/Eigenvalues> 
#include<fstream>
#include <iomanip>
#include<fastenvelope/FastEnvelope.h>// for is_triangle_degenerated()
namespace fastEnvelope {

	void obb::obb_init(const std::vector<std::vector<Vector3>>& envelope_vertices, const int p_face[8][3], const int c_face[6][3],
		const std::array<std::vector<int>, 8>& p_facepoint, const std::array<std::array<int, 4>, 6>& c_facepoint) {
		std::vector<Vector3> points;
		Trans.resize(envelope_vertices.size());
		invTrans.resize(envelope_vertices.size());
		int dege;
		Vector3 normal;

		for (int i = 0; i < envelope_vertices.size(); i++) {
			if (envelope_vertices[i].size() == polyhedron_point_number1) {
				Trans[i].resize(polyhedron_face_number1);
				invTrans[i].resize(polyhedron_face_number1);
				for (int j = 0; j < polyhedron_face_number1; j++) {
					//triangle to get normal: {envelope_vertices[i][p_face[j][0]],envelope_vertices[i][p_face[j][1]],envelope_vertices[i][p_face[j][2]]}
					dege = FastEnvelope::is_triangle_degenerated(envelope_vertices[i][p_face[j][0]], envelope_vertices[i][p_face[j][1]], envelope_vertices[i][p_face[j][2]]);
					if (dege == FastEnvelope::DEGENERATED_SEGMENT || dege == FastEnvelope::DEGENERATED_POINT) {
						std::cout << "need to fix here, face degeneration" << std::endl;
						exit(0);
					}
					
					if (dege == FastEnvelope::NERLY_DEGENERATED) {
						normal = FastEnvelope::accurate_normal_vector(envelope_vertices[i][p_face[j][0]], envelope_vertices[i][p_face[j][1]],
							envelope_vertices[i][p_face[j][0]], envelope_vertices[i][p_face[j][2]]);
					}
					if (dege == FastEnvelope::NOT_DEGENERATED) {
						normal = (envelope_vertices[i][p_face[j][0]] - envelope_vertices[i][p_face[j][1]]).cross(envelope_vertices[i][p_face[j][0]] - envelope_vertices[i][p_face[j][2]]).normalized();
					}
					points.clear();
					points.resize(p_facepoint[j].size());
					for (int k = 0; k < p_facepoint[j].size(); k++) {
						points[k] = envelope_vertices[i][p_facepoint[j][k]];
					}


					build_obb(points, normal, OBB_OFFSET, Trans[i][j], invTrans[i][j]);
				}

			}
			if (envelope_vertices[i].size() == polyhedron_point_number2) {
				Trans[i].resize(polyhedron_face_number2);
				invTrans[i].resize(polyhedron_face_number2);
				for (int j = 0; j < polyhedron_face_number2; j++) {
					//triangle to get normal: {envelope_vertices[i][p_face[j][0]],envelope_vertices[i][p_face[j][1]],envelope_vertices[i][p_face[j][2]]}
					dege = FastEnvelope::is_triangle_degenerated(envelope_vertices[i][p_face[j][0]], envelope_vertices[i][p_face[j][1]], envelope_vertices[i][p_face[j][2]]);
					if (dege == FastEnvelope::DEGENERATED_SEGMENT || dege == FastEnvelope::DEGENERATED_POINT) {
						std::cout << "need to fix here, face degeneration" << std::endl;
						exit(0);
					}

					if (dege == FastEnvelope::NERLY_DEGENERATED) {
						normal = FastEnvelope::accurate_normal_vector(envelope_vertices[i][p_face[j][0]], envelope_vertices[i][p_face[j][1]],
							envelope_vertices[i][p_face[j][0]], envelope_vertices[i][p_face[j][2]]);
					}
					if (dege == FastEnvelope::NOT_DEGENERATED) {
						normal = (envelope_vertices[i][p_face[j][0]] - envelope_vertices[i][p_face[j][1]]).cross(envelope_vertices[i][p_face[j][0]] - envelope_vertices[i][p_face[j][2]]).normalized();
					}
					points.clear();
					points.resize(p_facepoint[j].size());
					for (int k = 0; k < p_facepoint[j].size(); k++) {
						points[k] = envelope_vertices[i][p_facepoint[j][k]];
					}


					build_obb(points, normal, OBB_OFFSET, Trans[i][j], invTrans[i][j]);
				}
			}
		}
	}
	bool obb::intersected(const int prismid1, const int faceid1, const int prismid2, const int faceid2)const {
		return obb_intersection(Trans[prismid1][faceid1], invTrans[prismid1][faceid1], Trans[prismid2][faceid2], invTrans[prismid2][faceid2]);
	}
	void obb::build_obb(const std::vector<Vector3>& points, const Vector3& normal,const Scalar offset,
		Eigen::Matrix4d &M, Eigen::Matrix4d &invM)
	{
		Eigen::MatrixXd  PS(3, points.size());
		Vector3 cent;
		cent[0] = 0; cent[1] = 0; cent[2] = 0;
		for (int i = 0; i < points.size(); i++) {
			PS(0, i) = points[i][0];
			PS(1, i) = points[i][1];
			PS(2, i) = points[i][2];
			cent[0] += points[i][0];
			cent[1] += points[i][1];
			cent[2] += points[i][2];
		}
		cent[0] = cent[0] / points.size();
		cent[1] = cent[1] / points.size();
		cent[2] = cent[2] / points.size();
		Eigen::MatrixXd PCA= PS * PS.transpose();
		Eigen::EigenSolver<Eigen::MatrixXd> eg(PCA);
		Eigen::VectorXcd ev= eg.eigenvalues();
		
		int maxid = 0;

		for (int i = 1; i < 3; i++) {

			if (ev[maxid].real() <= ev(i).real()) {
				maxid = i;
			}
		}
		assert(ev(maxid).real() > 0);
		Eigen::VectorXcd v = eg.eigenvectors().col(maxid);
		Vector3 y;
		y[0] = v[0].real(); y[1] = v[1].real(); y[2] = v[2].real();
		y = y.normalized();
		assert(normal.dot(y) < SCALAR_ZERO);
		Vector3 z = normal.cross(y);
		Matrix3 R;
		R.col(0) = normal;
		R.col(1) = y;
		R.col(2) = z;

		Matrix3 invR = R.transpose();//inv(R)=trans(R)
		Vector3 T0 = -1 * R*cent, T1;//R*cent+T0=[0,0,0]
		Eigen::Matrix4d Trans0, invTrans0, Trans1, invTrans1, Scalling, invScalling;
		
		Trans0 << R, T0,
			0, 0, 0, 1;
		invTrans0 << invR, cent,
			0, 0, 0, 1;
		Eigen::MatrixXd  PSG(4, points.size()),prog(4, points.size());
		PSG << PS,
			Eigen::MatrixXd::Ones(1, points.size());
		prog = Trans0 * PSG;//contains the projection of the points in new axis
		Vector3 min, max, mid,corner;
		for (int i = 0; i < 3; i++) {
			min[i] = prog.row(i).minCoeff();
			max[i] = prog.row(i).maxCoeff();
		}
		mid = (min + max) / 2;
		corner = max - mid + Vector3(offset, offset, offset);
		Eigen::MatrixXd procg(4, 1), centng(4, 1), centn(3, 1);
		procg << mid,
			1;
		centng = invTrans0 * procg;
		centn(0) = centng(0); centn(1) = centng(1); centn(2) = centng(2);
		T1 = -1 * R*centn;
		Trans1 << R, T1,
			0, 0, 0, 1;
		invTrans1<< invR, centn,
			0, 0, 0, 1;
		Scalling = Eigen::Matrix4d::Zero(4, 4);
		Scalling(0, 0) = corner(0);
		Scalling(1, 1) = corner(1);
		Scalling(2, 2) = corner(2);
		Scalling(3, 3) = 1;
		invScalling = Scalling.inverse();
		M = invScalling * Trans1;
		invM = invTrans1 * Scalling;
	}
	
	bool obb::obb_intersection(const Eigen::Matrix4d &M1, const Eigen::Matrix4d &invM1, const Eigen::Matrix4d &M2, const Eigen::Matrix4d &invM2) {
		
		Eigen::MatrixXd ub(8, 4);
		
		ub << -1, -1, -1, 1,
			1, -1, -1, 1,
			1, 1, -1, 1,
			-1, 1, -1, 1,
			-1, 1, 1, 1,
			-1, -1, 1, 1,
			1, -1, 1, 1,
			1, 1, 1, 1;
		Eigen::MatrixXd b21(4, 8), b12(4, 8);
		
		Scalar minx, miny, minz, maxx, maxy, maxz;

		b12 = M2 * invM1*ub.transpose();

		minx = b12.row(0).minCoeff();
		if (minx >= 1) return false;
		miny= b12.row(1).minCoeff();
		if (miny >= 1) return false;
		minz = b12.row(2).minCoeff();
		if (minz >= 1) return false;

		maxx= b12.row(0).maxCoeff();
		if (maxx <= -1) return false;
		maxy = b12.row(1).maxCoeff();
		if (maxy <= -1) return false;
		maxz = b12.row(2).maxCoeff();
		if (maxz <= -1) return false;

		b21= M1 * invM2*ub.transpose();

		minx = b21.row(0).minCoeff();
		if (minx >= 1) return false;
		miny = b21.row(1).minCoeff();
		if (miny >= 1) return false;
		minz = b21.row(2).minCoeff();
		if (minz >= 1) return false;

		maxx = b21.row(0).maxCoeff();
		if (maxx <= -1) return false;
		maxy = b21.row(1).maxCoeff();
		if (maxy <= -1) return false;
		maxz = b21.row(2).maxCoeff();
		if (maxz <= -1) return false;


		return true;
		
	}
#include <time.h>
	void obb::test(){
		std::vector<Vector3> p, p1;
		Vector3 dis;
		srand(int(time(0)));
		dis = Vector3(rand() % 100, rand() % 100, rand() % 100) / 50;
		p.resize(20);
		for (int i = 0; i < p.size(); i++) {
			p[i] = Vector3(rand() % 100, rand() % 100, rand() % 100) / 100 + dis;
		}
		Eigen::MatrixXd  PS(3, p.size());
		

		for (int i = 0; i < p.size(); i++) {
			PS(0, i) = p[i][0];
			PS(1, i) = p[i][1];
			PS(2, i) = p[i][2];
			
		}
		
		Eigen::MatrixXd PCA = PS * PS.transpose();
		Eigen::EigenSolver<Eigen::MatrixXd> eg(PCA);
		Eigen::VectorXcd ev = eg.eigenvalues();

		int minid = 0;

		for (int i = 1; i < 3; i++) {

			if (ev[minid].real() >= ev(i).real()) {//find the minimal eigenvalue
				minid = i;
			}
		}

		Eigen::VectorXcd v = eg.eigenvectors().col(minid);
		Vector3 normal;
		normal[0] = v[0].real(); normal[1] = v[1].real(); normal[2] = v[2].real();
		Scalar offset = 0.1;
		Eigen::Matrix4d M, invM;
		build_obb(p, normal, offset, M, invM);
		std::ofstream fout;
		fout.open("D:\\vs\\fast_envelope\\obb\\points.txt");
		
		for (int i = 0; i < p.size(); i++) {

			fout << std::setprecision(17) << p[i][0] << " " << p[i][1] << " " << p[i][2] << std::endl;

		}
		fout.close();
		fout.open("D:\\vs\\fast_envelope\\obb\\matrixes.txt");

		for (int i = 0; i < 4; i++) {

			fout << std::setprecision(17) << M(i, 0) << " " << M(i, 1) << " " << M(i, 2) << " " << M(i, 3) << std::endl;

		}
		for (int i = 0; i < 4; i++) {

			fout << std::setprecision(17) << invM(i, 0) << " " << invM(i, 1) << " " << invM(i, 2) << " " << invM(i, 3) << std::endl;

		}

		fout.close();




		//////////////////////////
	
		//srand(int(time(0)));
		dis = Vector3(rand() % 100, rand() % 100, rand() % 100) / 100*3;
		p.clear();
		p.resize(20);
		for (int i = 0; i < p.size(); i++) {
			p[i] = Vector3(rand() % 100, rand() % 100, rand() % 100) / 100 + dis;
		}
		Eigen::MatrixXd  PS1(3, p.size());


		for (int i = 0; i < p.size(); i++) {
			PS1(0, i) = p[i][0];
			PS1(1, i) = p[i][1];
			PS1(2, i) = p[i][2];

		}

		Eigen::MatrixXd PCA1 = PS1 * PS1.transpose();
		Eigen::EigenSolver<Eigen::MatrixXd> eg1(PCA1);
		Eigen::VectorXcd ev1 = eg1.eigenvalues();

		minid = 0;

		for (int i = 1; i < 3; i++) {

			if (ev1[minid].real() >= ev1(i).real()) {//find the minimal eigenvalue
				minid = i;
			}
		}

		Eigen::VectorXcd v1 = eg1.eigenvectors().col(minid);
		
		normal[0] = v1[0].real(); normal[1] = v1[1].real(); normal[2] = v1[2].real();
		
		Eigen::Matrix4d M1, invM1;
		build_obb(p, normal, offset, M1, invM1);

		fout.open("D:\\vs\\fast_envelope\\obb\\points1.txt");

		for (int i = 0; i < p.size(); i++) {

			fout << std::setprecision(17) << p[i][0] << " " << p[i][1] << " " << p[i][2] << std::endl;

		}
		fout.close();
		fout.open("D:\\vs\\fast_envelope\\obb\\matrixes1.txt");

		for (int i = 0; i < 4; i++) {

			fout << std::setprecision(17) << M1(i, 0) << " " << M1(i, 1) << " " << M1(i, 2) << " " << M1(i, 3) << std::endl;

		}
		for (int i = 0; i < 4; i++) {

			fout << std::setprecision(17) << invM1(i, 0) << " " << invM1(i, 1) << " " << invM1(i, 2) << " " << invM1(i, 3) << std::endl;

		}

		fout.close();

		std::cout << "is intersected? " << obb_intersection(M, invM, M1, invM1) << std::endl;
	}
}