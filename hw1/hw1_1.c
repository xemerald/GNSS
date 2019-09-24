#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <constant.h>
#include <formulas.h>

int main( int argc, char **argv )
{
	double orb_height = 0.0;
	double res = 0.0;

	if ( argc < 2 ) {
		printf("Usage: %s Orbis_Height\n", argv[0]);
		return -1;
	}

	orb_height = atof(argv[1]) + EARTH_RAD;
/* Derive the circular orbit period */
	res = get_cir_orbit_period( orb_height );
	printf("Period: %lf sec (~ %lf min).\n", res, res/60.0);
/* Derive the circular orbit velocity */
	res = get_cir_orbit_velocity( orb_height );
	printf("Velocity: %lf m/s (~ %lf km/s).\n", res, res/1000.0);

	return 0;
}
