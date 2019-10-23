#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* PGPLOT C library header */
#include <cpgplot.h>
#include <constant.h>
#include <formulas.h>
#include <matrix.h>

int main( int argc, char **argv )
{
	int    i;
	double origin_lat, origin_lon;
	double input_y[5] = { 0.0 };
	double input_c[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
	double input_x[5] = { 2015.0, 2016.0, 2017.0, 2018.0, 2019.0 };
	double output[2] = { 0.0 };
	float y[5] = { 0.0 };
	float year[5] = { 2015.0, 2016.0, 2017.0, 2018.0, 2019.0 };
	float max_y = 0.0;
	char testc[120];

	MATRIX *d = NULL;
	MATRIX *g = NULL;
	MATRIX *m = NULL;

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
	}

	d = matrix_new( 5, 1 );
	g = matrix_new( 5, 2 );

	matrix_assign_col( d, input_y, 1, 5 );
	matrix_assign_col( g, input_c, 1, 5 );
	matrix_assign_col( g, input_x, 2, 5 );

	m = matrix_div( d, g );
	matrix_free( d );
	matrix_free( g );
	matrix_extract_seq( m, output, 2 );
	matrix_free( m );
	printf("%lf %lf\n", output[0], output[1]);
/*  */

/* Just plot the result by PGPLOT library */
	if ( !cpgopen("position_with_time/cps") ) return -1;
	cpgenv( 2014.0, 2020.0, 0.0, max_y*1.01, 0, 1);
	cpglab("Year", "Position(cm)", "Position Change over Years in Lon. (st. S01R)");
	cpgsci(2);
	cpgpt(5, year, y, 3);
	cpgsci(1);
	cpgmove(2014.0, output[0] + output[1] * 2014.0);
	cpgdraw(2020.0, output[0] + output[1] * 2020.0);
	sprintf(testc, "Y = %.3f + %.3f * X", output[0], output[1]);
	cpgtext(2014.5, 10.0, testc);
	sprintf(testc, "Velocity: %.3f cm/yr", output[1]);
	cpgtext(2014.5, 9.5, testc);
	cpgend();


	return 0;
}
