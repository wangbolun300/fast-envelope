#pragma once
#include<math.h>
int get_projection_plane(
    double ov1x, double ov1y, double ov1z,
    double ov2x, double ov2y, double ov2z,
    double ov3x, double ov3y, double ov3z);

void triangle_normal_exact(
	double ov1x, double ov1y, double ov1z,
	double ov2x, double ov2y, double ov2z,
	double ov3x, double ov3y, double ov3z,
	double &nvxc, double &nvyc, double &nvzc);

int dot_product_sign(double vx, double vy, double vz,
	double ux, double uy, double uz);

bool is_tpp_in_triangle(double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
	double w1x, double w1y, double w1z, double w2x, double w2y, double w2z, double w3x, double w3y, double w3z,
	double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double u3x, double u3y, double u3z,
	double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double q3x, double q3y, double q3z);


