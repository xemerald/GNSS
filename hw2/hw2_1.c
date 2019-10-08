#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <constant.h>
#include <formulas.h>

int main( int argc, char **argv )
{
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	GEOCORD res;

	if ( argc < 4 ) {
		printf("Usage: %s X Y Z\n", argv[0]);
		return -1;
	}

	x = atof(argv[1]);
	y = atof(argv[2]);
	z = atof(argv[3]);
/* Derive the circular orbit period */
	res = xyz2geocd( x, y, z );
	printf("Geo Position: Longitude: %lf, Latitude: %lf\n", res.longitude, res.latitude);

	return 0;
}
