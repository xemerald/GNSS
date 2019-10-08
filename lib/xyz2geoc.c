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
