#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <constant.h>
#include <formulas.h>
#include <matrix.h>

/**/
static int arg_parse( int, char ** );
static double *fill_value( double *, double, double, double );

int main( int argc, char **argv )
{
	int    i, axis;
	double Sxx, Sxy, Sxz, Syy, Syz, Szz;
	double ang;
	double input_array[3];
	double output_array[9];

	MATRIX *d  = NULL;
	MATRIX *s  = NULL;
	MATRIX *dt = NULL;
	MATRIX *ds = NULL;
	MATRIX *dsdt = NULL;

/**/
	if ( arg_parse( argc, argv ) ) return -1;

	s = matrix_new( 3, 3 );
	d = matrix_new( 3, 3 );

	Sxx  = atof(argv[1]);
	Sxy  = atof(argv[2]);
	Sxz  = atof(argv[3]);
	Syy  = atof(argv[4]);
	Syz  = atof(argv[5]);
	Szz  = atof(argv[6]);
	axis = atoi(argv[7]);
	ang  = atof(argv[8]) * PI / 180.0;

	fill_value( input_array, Sxx, Sxy, Sxz );
	matrix_assign_row( s, input_array, 1, 3 );
	fill_value( input_array, Sxy, Syy, Syz );
	matrix_assign_row( s, input_array, 2, 3 );
	fill_value( input_array, Sxz, Syz, Szz );
	matrix_assign_row( s, input_array, 3, 3 );

	switch ( axis ) {
		case 0:
			fill_value( input_array, Sxx, Sxy, Sxz );
			matrix_assign_row( d, input_array, 1, 3 );
			fill_value( input_array, Sxy, Syy, Syz );
			matrix_assign_row( d, input_array, 2, 3 );
			fill_value( input_array, Sxz, Syz, Szz );
			matrix_assign_row( d, input_array, 3, 3 );
			break;
		case 1:
			fill_value( input_array, cos(ang), 0.0, cos(0.5*PI + ang) );
			matrix_assign_row( d, input_array, 1, 3 );
			fill_value( input_array, 0.0, 1.0, 0.0 );
			matrix_assign_row( d, input_array, 2, 3 );
			fill_value( input_array, cos(0.5*PI - ang), 0.0, cos(ang) );
			matrix_assign_row( d, input_array, 3, 3 );
			break;
		case 2:
			fill_value( input_array, cos(ang), cos(0.5*PI + ang), 0.0 );
			matrix_assign_row( d, input_array, 1, 3 );
			fill_value( input_array, cos(0.5*PI - ang), cos(ang), 0.0 );
			matrix_assign_row( d, input_array, 2, 3 );
			fill_value( input_array, 0.0, 0.0, 1.0 );
			matrix_assign_row( d, input_array, 3, 3 );
			break;
		default:
			break;
	}

	dt = matrix_transpose( d );
	ds = matrix_mul( d, s );
	dsdt = matrix_mul( ds, dt );

	matrix_extract_seq( dsdt, output_array, 9 );

	for ( i=0; i<9; i++ ) printf("%e ", output_array[i]);
	printf("\n");

	matrix_free( d );
	matrix_free( s );
	matrix_free( ds );
	matrix_free( dt );
	matrix_free( dsdt );

	return 0;
}

/**/
static int arg_parse( int argc, char **argv ) {

/**/
	if ( argc < 9 ) {
		fprintf(stderr, "Usage: %s Sxx Sxy Sxz Syy Syz Szz Rotation_axis Rotation_angle\n", argv[0]);
		return -1;
	}

	return 0;
}

static double *fill_value( double *dest, double val1, double val2, double val3 ) {

	dest[0] = val1;
	dest[1] = val2;
	dest[2] = val3;

	return dest;
}
