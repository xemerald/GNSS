#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* PGPLOT C library header */
#include <cpgplot.h>
#include <constant.h>
#include <formulas.h>

int main( int argc, char **argv )
{
	int    i;
	double origin_lat, origin_lon;
	double input_y[5] = { 0.0 };
	double input_x[5] = { 2015.0, 2016.0, 2017.0, 2018.0, 2019.0 };
	float y[5] = { 0.0 };
	float year[5] = { 2015.0, 2016.0, 2017.0, 2018.0, 2019.0 };
	float max_y = 0.0;
	char testc[120];


	if ( argc < 6 ) {
		printf("Usage: %s P1 P2 P3 P4 P5\n", argv[0]);
		return -1;
	}

	for ( i=0; i < 5; i++ )
		input_y[i] = atof(argv[i+1]);

	origin_lon = input_y[0];
	origin_lat = 23.655273;

	for ( i=0; i < 5; i++ ) {
		//input_y[i] = coor2distf( origin_lat, origin_lon, input_y[i], origin_lon );
		input_y[i] = coor2distf( origin_lat, origin_lon, origin_lat, input_y[i] );
		y[i] = input_y[i];
		if ( y[i] > max_y ) max_y = y[i];
		printf("%f\n", y[i]);
	}

	SLinearRegression(input_x, input_y, 5, &origin_lat, &origin_lon);
	printf("%f %f\n", origin_lat, origin_lon);
/*  */

/* Just plot the result by PGPLOT library */
	if ( !cpgopen("position_with_time/cps") ) return -1;
	cpgenv( 2014.0, 2020.0, 0.0, max_y*1.01, 0, 1);
	cpglab("Year", "Position(cm)", "Position Change over Years in Lon. (st. S01R)");
	cpgsci(2);
	cpgpt(5, year, y, 3);
	cpgsci(1);
	cpgmove(2014.0, origin_lat + origin_lon * 2014.0);
	cpgdraw(2020.0, origin_lat + origin_lon * 2020.0);
	sprintf(testc, "Y = %.3f + %.3f * X", origin_lat, origin_lon);
	cpgtext(2014.5, 10.0, testc);
	sprintf(testc, "Velocity: %.3f cm/yr", origin_lon);
	cpgtext(2014.5, 9.5, testc);
	cpgend();


	return 0;
}
