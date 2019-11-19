#pragma once

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
