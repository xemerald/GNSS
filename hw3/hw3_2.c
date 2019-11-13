#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <constant.h>
#include <formulas.h>
#include <matrix.h>

/**/
static int arg_parse( int, char ** );

int main( int argc, char **argv )
{
	int     i;
	double alat, alon, blat, blon;
	double avelx, avely;
	double delta_x, delta_y;
	double input_array[4];

	MATRIX *d  = NULL;
	MATRIX *g  = NULL;
	MATRIX *ig  = NULL;
	MATRIX *m  = NULL;

/**/
	if ( arg_parse( argc, argv ) ) return -1;

	g = matrix_new( 4, 4 );
	d = matrix_new( 4, 1 );

	alat  = atof(argv[1]);
	alon  = atof(argv[2]);
	avelx = atof(argv[7]);
	avely = atof(argv[8]);
	matrix_prefill_array( input_array, 4,
		avelx - atof(argv[3]), avely - atof(argv[4]), atof(argv[11]) - avelx, atof(argv[12]) - avelx );
	matrix_assign_col( d, input_array, 1, 4 );

	blat = atof(argv[5]);
	blon = atof(argv[6]);
	delta_x = coor2distf( blat, alon, blat, blon ) * 10.0 * ( blon > alon ? 1.0 : -1.0 ); /* unit: mm */
	delta_y = coor2distf( alat, blon, blat, blon ) * 10.0 * ( blat > alat ? 1.0 : -1.0 ); /* unit: mm */
	matrix_prefill_array( input_array, 4, delta_x, delta_y, 0.0, delta_y );
	matrix_assign_row( g, input_array, 1, 4 );
	matrix_prefill_array( input_array, 4, 0.0, delta_x, delta_y, -delta_x );
	matrix_assign_row( g, input_array, 2, 4 );

	alat = atof(argv[9]);
	alon = atof(argv[10]);
	delta_x = coor2distf( blat, alon, blat, blon ) * 10.0 * ( blon < alon ? 1.0 : -1.0 ); /* unit: mm */
	delta_y = coor2distf( alat, blon, blat, blon ) * 10.0 * ( blat < alat ? 1.0 : -1.0 ); /* unit: mm */
	matrix_prefill_array( input_array, 4, delta_x, delta_y, 0.0, delta_y );
	matrix_assign_row( g, input_array, 3, 4 );
	matrix_prefill_array( input_array, 4, 0.0, delta_x, delta_y, -delta_x );
	matrix_assign_row( g, input_array, 4, 4 );

	ig = matrix_inverse( g );
	m  = matrix_mul( ig, d );

	matrix_extract_seq( m, input_array, 4 );

	for ( i=0; i<4; i++ ) printf("%.3e strain\n", input_array[i]);

	matrix_free( d );
	matrix_free( ig );
	matrix_free( g );
	matrix_free( m );

	return 0;
}

/**/
static int arg_parse( int argc, char **argv ) {

/**/
	if ( argc < 13 ) {
		fprintf(stderr, "Usage: %s lat_sta1 long_sta1 vx_sta1 vy_sta1 lat_sta2 long_sta2 vx_sta2 vy_sta2 lat_sta3 long_sta3 vx_sta3 vy_sta3\n", argv[0]);
		return -1;
	}

	return 0;
}
