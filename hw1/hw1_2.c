#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* PGPLOT C library header */
#include <cpgplot.h>
#include <constant.h>
#include <formulas.h>

int main( int argc, char **argv )
{
	int i, max_i;
	double orb_height = 0.0;
	double radius = 0.0;
	double freq = 0.0;

	double delta = 0.0;
	double deg_inc = 0.0, theta = 0.0;
	float *x = NULL, *y = NULL;
	float max_y = 0.0, min_y = 0.0;

	if ( argc < 3 ) {
		printf("Usage: %s Orbis_Height Signal_Frequency\n", argv[0]);
		return -1;
	}

/* Setting the parameters */
	orb_height = atof(argv[1]);
	radius     = orb_height + EARTH_RAD;
	freq       = atof(argv[2]);
/* Get the period of this orbit & derive the degree incremental in one time (one sencond) step */
	delta      = get_cir_orbit_period(radius);
	deg_inc    = 360.0 / delta;
	delta      = radius * radius;
/* Derive the circular angle of observable arc */
	theta      = asin((delta - EARTH_RAD * EARTH_RAD) / delta) * 180.0 / PI;
/* Allocate the memory space for results */
	max_i      = 4096;
	x          = calloc(max_i, sizeof(float));
	y          = calloc(max_i, sizeof(float));

/* Compute the data series along with time, I assume that it should be equal angular-velocity motion
   so the incremental degree within the specified time interval should be the same all the way. */
	for ( delta = theta, i = 0; delta >= -theta; delta -= deg_inc, i++ ) {
	/* Simply use the formula to derive the frequency shift at this position. And the input circular angle(delta) is
	   the angle between the satellite and the zenith of receiver. */
		y[i] = get_freq_shift_cir_angle(delta, orb_height, freq) - freq;
	/* Each step should be one second here. */
		x[i] = i * 1.0;
	/* Try to find out the maximum and minimum value of frequency shift. It would be the y-axis limits when plotting. */
		if ( y[i] > max_y ) max_y = y[i];
		if ( y[i] < min_y ) min_y = y[i];

	/* If the number of sample is larger than previous setting, it would be re-allocated here */
		if( i >= max_i ) {
			max_i += 1024;
			x = realloc(x, max_i * sizeof(float));
			y = realloc(y, max_i * sizeof(float));
		}
	}

/* Just plot the result by PGPLOT library */
	if ( !cpgopen("freq_shift_time/cps") ) return -1;
	cpgenv( 0.0, 2.0*theta/deg_inc, min_y*1.01, max_y*1.01, 0, 1);
	cpglab("Time(sec)", "Frequency shift(Hz)", "Freq. shift along time");
	cpgline( i, x, y );
	cpgend();

/* Compute the data series along with observing angle; therefore, the incremental unit should be observing elevation angle
   in degree here. And what we can see should always above the horizon(0 to 180 degree). */
	for ( delta = 0.0, i = 0; delta <= 180.0; delta += 1.0, i++ ) {
	/* Here would be much simpler than previous condition that the input angle(delta) is just the elevation angle between
	   the horizon and the satellite. */
		y[i] = get_freq_shift_obs_angle(delta, orb_height, freq) - freq;
	/* Each step should be one degree. */
		x[i] = i * 1.0;
		if ( y[i] > max_y ) max_y = y[i];
		if ( y[i] < min_y ) min_y = y[i];

		if( i >= max_i ) {
			max_i += 1024;
			x = realloc(x, max_i * sizeof(float));
			y = realloc(y, max_i * sizeof(float));
		}
	}

	if ( !cpgopen("freq_shift_angle/cps") ) return -1;
	cpgenv( 0.0, 180.0, min_y*1.01, max_y*1.01, 0, 1);
	cpglab("Obs. angle(deg)", "Frequency shift(Hz)", "Freq. shift along obs. angle");
	cpgline( i, x, y );
	cpgend();

	return 0;
}
