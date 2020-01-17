#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "../include/constant.h"
#include "../include/faultgeom.h"
#include "../include/matrix.h"

static double add_value( const double );

/*
	Move the fault so that the coordinates of the midpoint
	refer to the fault bottom.
*/
FAULT_GEOM transfault_bottom( FAULT_GEOM input )
{
	const double deg2rad = DEG2RAD;
	const double diprad  = input.dip * deg2rad;

	double dipdir = (input.strike + 90.0) * deg2rad;
	double offset = fabs(input.width * cos(diprad));

	input.depth += fabs(input.width * sin(diprad));
	input.dip   += 180.0;
	input.east  += offset * sin(dipdir);
	input.north += offset * cos(dipdir);

	return input;
}

static double _value = 0.0;

/*
	This functions discretizes a fault model into i*j distinct patches.

	INPUTS:
    	input = 1x7 vector defined as follows
    	i = number of patches along fault length
    	j = number of patches along fault width

	OUTPUTS:
    	patches = mx7 matrix of patch models, where m=i*j
*/
int fault2patch( const FAULT_GEOM input, const int i, const int j, FAULT_GEOM **patches )
{
	int _i, _j;

	const double deg2rad   = DEG2RAD;
	const double diprad    = input.dip * deg2rad;
	const double strikerad = -input.strike * deg2rad;
	const double sin_dip   = sin(diprad);
	const double cos_dip   = cos(diprad);
	const double iwidth    = input.length / (double)i;
	const double jwidth    = input.width  / (double)j;
	const int    n         = i * j;

	double tmp = 0.0;
	double *fill = (double *)malloc((3*n > 9 ? 3*n : 9)*sizeof(double));
	double *fillptr = NULL;
	MATRIX *a_mat, *b_mat;
	MATRIX *mp;

	FAULT_GEOM *output = (FAULT_GEOM *)calloc(sizeof(FAULT_GEOM), n);

	a_mat = matrix_new( n, 3 );

/* */
	for ( _j = 0, tmp = jwidth; _j < j; _j++, tmp += jwidth )
		fill[_j] = cos_dip * (tmp - input.width);
/* */
	for ( _i = 1; _i < i; _i++ )
		memcpy(fill + j*_i, fill, j*sizeof(double));
/* */
	matrix_assign_col( a_mat, fill, 1, n );

/* */
	for ( _j = 0, tmp = jwidth * (j - 1); _j < j; _j++, tmp -= jwidth )
		fill[_j] =  input.depth - sin_dip * tmp;
/* */
	for ( _i = 1; _i < i; _i++ )
		memcpy(fill + j*_i, fill, j*sizeof(double));
/* */
	matrix_assign_col( a_mat, fill, 3, n );

/* */
	fillptr = fill;
	for ( _i = 0, tmp = iwidth; _i < i; _i++, tmp += iwidth ) {
		*fillptr = tmp - 0.5 * (input.length + iwidth);
		for ( _j = 1; _j < j; _j++ ) *fillptr = *fillptr++;
	}
/* */
	matrix_assign_col( a_mat, fill, 2, n );

/* */
	b_mat = matrix_new( 3, 3 );
	matrix_prefill_array( fill, 9, cos(strikerad), sin(strikerad), 0.0, -sin(strikerad), cos(strikerad), 0.0, 0.0, 0.0, 1.0 );
	matrix_assign_seq( b_mat, fill, 9 );

	mp = matrix_mul( a_mat, b_mat );
	matrix_free( a_mat );
	matrix_free( b_mat );

	_value = input.east;
	matrix_apply_col( mp, add_value, 1 );
	_value = input.north;
	matrix_apply_col( mp, add_value, 2 );

	matrix_extract_seq(mp, fill, 3*n);
	matrix_free(mp);
/* */
	fillptr = fill;
	for( _i = 0; _i < n; _i++ ) {
		output[_i].length = iwidth;
		output[_i].width  = jwidth;
		output[_i].dip    = input.dip;
		output[_i].strike = input.strike;
		output[_i].east   = *fillptr++;
		output[_i].north  = *fillptr++;
		output[_i].depth  = *fillptr++;
	}

	free(fill);
	*patches = output;

	return n;
}

static double add_value( const double input )
{
	return input + _value;
}
