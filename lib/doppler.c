#include <math.h>
#include "../include/formulas.h"
#include "../include/constant.h"

/* The basic doppler effect formula, it can output the shifted frequency. */
static double simple_doppler_func( double wa_speed, double re_speed, double so_speed, double freq )
{
	return ((wa_speed + re_speed) / (wa_speed - so_speed) * freq);
}

/* Base on the doppler effect compute the shifted frequency of the satellite signal by
   the receiver observing elevation angle. */
double get_freq_shift_obs_angle( double obs_angle, double orb_height, double freq )
{
/* Turn it into orbit radius */
	orb_height += EARTH_RAD;
	double velocity = get_cir_orbit_velocity( orb_height );
/* Plus 90 deg. & convert degree to radian */
	obs_angle = (obs_angle + 90.0) * PI / 180.0;
/* We just assume that the satellite orbit should always across the zenith of receiver and the receiver is
   fixed(no related motion). And once the observing elevation angle growth the angle between two vector:
   satellite to earth center and satellite to receiver will decrease until over the zenith. The angle between
   the two vector is the key to derive the receiver-oriented velocity of satellite due to the this velocity is
   equal the the tangent velocity of satellite multiply by the sine of the angle. Since that we can derive the
   angle by the law of sines so the result is shown as follow... */
	velocity *= sin(obs_angle) * EARTH_RAD / orb_height;

/* Once we get the relative velocity between satellite and receiver, we simply use the doppler formula to get the
   shifted frequency. */
	return simple_doppler_func( LIGHT_SPEED, 0.0, velocity, freq );
}

/* Base on the doppler effect compute the shifted frequency of the satellite signal by
   the circular angle of the arc between the satellite and the zenith of receiver. */
double get_freq_shift_cir_angle( double cir_angle, double orb_height, double freq )
{
/* Turn it into orbit radius */
	orb_height += EARTH_RAD;
	double velocity = get_cir_orbit_velocity( orb_height );
/* Convert degree to radian */
	cir_angle *= PI / 180.0;
/* Again, we just assume that the satellite orbit should always across the zenith of receiver and the receiver is
   fixed(no related motion). But this time we will use the circular angle of the arc between the satellite and the zenith
   of receiver to derive the receiver-oriented velocity of satellite. However, before that we still need the angle between
   two vector: satellite to earth center and satellite to receiver. And the relationship between the circular angle and
   the angle between the two vector can be described by the law of sines and the law of cosines so the result is shown
   as follow... */
	velocity  *= sin(cir_angle) * EARTH_RAD / sqrt(EARTH_RAD*EARTH_RAD + orb_height*orb_height - 2.0*EARTH_RAD*orb_height*cos(cir_angle));

/* Once we get the relative velocity between satellite and receiver, we simply use the doppler formula to get the
   shifted frequency. */
	return simple_doppler_func( LIGHT_SPEED, 0.0, velocity, freq );
}
