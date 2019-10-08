#include <stdio.h>
#include <math.h>
#include "../include/formulas.h"
#include "../include/constant.h"

/*  */
GEOCORD xyz2geocd( double x, double y, double z )
{
	const double a = EARTH_LAXIS;
	const double b = a - a / EARTH_INVFLAT;
	GEOCORD ret    = { 0.0, atan2(y, x), sqrt(x*x + y*y) };
	double beta    = atan2(a*z, b*ret.elevation);

	x = sin(beta);
	x *= x*x;
	y = cos(beta);
	y *= y*y;
	beta = 180.0 / PI;

	ret.latitude   = atan2(z + (a*a/b - b)*x, ret.elevation - (a - b*b/a)*y) * beta;
	ret.longitude *= beta;

	return ret;
}

/***********************************************************************************
 * coor2distf() Transforms the coordinate(latitude & longitude) into distance(unit: cm) *
 ***********************************************************************************/
double coor2distf( const double elat, const double elon, const double slat, const double slon )
{
	const double avlat = (elat + slat)*0.5;

	double a = 1.840708 + avlat*(.0015269 + avlat*(-.00034 + avlat*(1.02337e-6)));
	double b = 1.843404 + avlat*(-6.93799e-5 + avlat*(8.79993e-6 + avlat*(-6.47527e-8)));

	a *= (slon - elon) * 60.0;
	b *= (slat - elat) * 60.0;

	return sqrt(a*a + b*b) * 100000.0;
}


int SLinearRegression(double *xi, double *yi, int rows, double *a, double *b)
{
	int m;
	double ux = 0.0, uy = 0.0;
	double Lxx = 0.0, Lxy = 0.0;

	if ( xi == NULL || yi == NULL || a == NULL || b == NULL || rows < 1 ) return -1;

	for ( m = 0; m < rows; m++ ) {
		ux += *(xi + m);
		uy += *(yi + m);
	}

	ux /= (double)rows;
	uy /= (double)rows;

	for ( m = 0; m < rows; m++ ) {
		Lxx += (*(xi + m) - ux)*(*(xi + m) - ux);
		Lxy += (*(xi + m) - ux)*(*(yi + m) - uy);
	}

	*b = Lxy / Lxx;
	*a = uy - *b * ux;

	return 0;
}
