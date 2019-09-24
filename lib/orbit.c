#include <math.h>
#include "../include/constant.h"

/* Simply use the formula: t^2 / a^3 = 4*pi / G*M */
double get_cir_orbit_period( double r )
{
	return sqrt((4.0 * PI_SQ) / EARTH_MASS_G * r * r * r);
}

/* Simply use the formula: v^2 = G*M / R */
double get_cir_orbit_velocity( double r )
{
	return sqrt(EARTH_MASS_G / r);
}
