#pragma once

typedef struct {
	double latitude;
	double longitude;
	double elevation;
} GEOCORD;

double get_cir_orbit_period( double );
double get_cir_orbit_velocity( double );

double get_freq_shift_obs_angle( double, double, double );
double get_freq_shift_cir_angle( double, double, double );

GEOCORD xyz2geocd( double, double, double );
double coor2distf( const double, const double, const double, const double );
int SLinearRegression(double *, double *, int, double *, double *);
