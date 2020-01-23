int get_projection_plane(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z)
{
    int r = triangle_normal_filtered(ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x, ov3y, ov3z);
    if (r >= 0)
        return r;

    return triangle_normal_exact(ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x, ov3y, ov3z);
}

void triangle_normal_exact(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z, double &nvxc, double &nvyc, double &nvzc)
{
    expansionObject o;
    double v3x[2];
    o.two_Diff(ov3x, ov2x, v3x);
    double v3y[2];
    o.two_Diff(ov3y, ov2y, v3y);
    double v3z[2];
    o.two_Diff(ov3z, ov2z, v3z);
    double v2x[2];
    o.two_Diff(ov2x, ov1x, v2x);
    double v2y[2];
    o.two_Diff(ov2y, ov1y, v2y);
    double v2z[2];
    o.two_Diff(ov2z, ov1z, v2z);
    double nvx1[8];
    o.Two_Two_Prod(v2y, v3z, nvx1);
    double nvx2[8];
    o.Two_Two_Prod(v2z, v3y, nvx2);
    double nvx[16];
    int nvx_len = o.Gen_Diff(8, nvx1, 8, nvx2, nvx);
    double nvy1[8];
    o.Two_Two_Prod(v3x, v2z, nvy1);
    double nvy2[8];
    o.Two_Two_Prod(v3z, v2x, nvy2);
    double nvy[16];
    int nvy_len = o.Gen_Diff(8, nvy1, 8, nvy2, nvy);
    double nvz1[8];
    o.Two_Two_Prod(v2x, v3y, nvz1);
    double nvz2[8];
    o.Two_Two_Prod(v2y, v3x, nvz2);
    double nvz[16];
    int nvz_len = o.Gen_Diff(8, nvz1, 8, nvz2, nvz);

    nvxc = nvx[nvx_len - 1];
    nvyc = nvy[nvy_len - 1];
    nvzc = nvz[nvz_len - 1];

    double l = sqrt(nvxc * nvxc + nvyc * nvyc + nvzc * nvzc);

    nvxc /= l;
    nvyc /= l;
    nvzc /= l;
}

void cross_product_normalized_exact(
	double ov1x, double ov1y, double ov1z,
	double ov2x, double ov2y, double ov2z,
	double pv1x, double pv1y, double pv1z,
	double pv2x, double pv2y, double pv2z,
	double &nvxc, double &nvyc, double &nvzc) {
	expansionObject o;
	double v3x[2];
	o.two_Diff(pv2x, pv1x, v3x);
	double v3y[2];
	o.two_Diff(pv2y, pv1y, v3y);
	double v3z[2];
	o.two_Diff(pv2z, pv1z, v3z);
	double v2x[2];
	o.two_Diff(ov2x, ov1x, v2x);
	double v2y[2];
	o.two_Diff(ov2y, ov1y, v2y);
	double v2z[2];
	o.two_Diff(ov2z, ov1z, v2z);
	double nvx1[8];
	o.Two_Two_Prod(v2y, v3z, nvx1);
	double nvx2[8];
	o.Two_Two_Prod(v2z, v3y, nvx2);
	double nvx[16];
	int nvx_len = o.Gen_Diff(8, nvx1, 8, nvx2, nvx);
	double nvy1[8];
	o.Two_Two_Prod(v3x, v2z, nvy1);
	double nvy2[8];
	o.Two_Two_Prod(v3z, v2x, nvy2);
	double nvy[16];
	int nvy_len = o.Gen_Diff(8, nvy1, 8, nvy2, nvy);
	double nvz1[8];
	o.Two_Two_Prod(v2x, v3y, nvz1);
	double nvz2[8];
	o.Two_Two_Prod(v2y, v3x, nvz2);
	double nvz[16];
	int nvz_len = o.Gen_Diff(8, nvz1, 8, nvz2, nvz);

	nvxc = nvx[nvx_len - 1];
	nvyc = nvy[nvy_len - 1];
	nvzc = nvz[nvz_len - 1];

	double l = sqrt(nvxc * nvxc + nvyc * nvyc + nvzc * nvzc);

	nvxc /= l;
	nvyc /= l;
	nvzc /= l;
}


int dot_product_sign(double vx, double vy, double vz,
	double ux, double uy, double uz) {
	expansionObject o;
	double x[2];
	o.Two_Prod(vx, ux, x);
	double y[2];
	o.Two_Prod(vy, uy, y);
	double z[2];
	o.Two_Prod(vz, uz, z);
	double xsumy[4];
	o.Two_Two_Sum(x[1], x[0], y[1], y[0], xsumy[3], xsumy[2], xsumy[1], xsumy[0]);
	double sum[6];
	int sl = o.Gen_Sum(4, xsumy, 2, z, sum);
	double s = sum[sl - 1];

	if (s > 0) return 1;
	if (s < 0) return -1;
	return 0;
}
