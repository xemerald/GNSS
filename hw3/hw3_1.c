#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* PGPLOT C library header */
#include <cpgplot.h>
#include <constant.h>
#include <formulas.h>
#include <matrix.h>

#define TAU          3.0f
#define NUM_PARAM    6
#define JIASHIAN_EQ  2010.1726f
#define MEINONG_EQ   2016.1014f

static double seasonal_term_sin( const double );
static double seasonal_term_cos( const double );
static double heavy_side( const double );
static double post_seismic_def( const double );
static double heavy_side_eq_jiashian( const double );
static double post_seismic_def_eq_jiashian( double );
static double final_function( const double, const double, const double, const double, const double, const double, const double );

int main( int argc, char **argv )
{
	int    i, j, count = 0;
	double pivot[3] = { 0.0 };
	double min_y, max_y;
	double *input_array  = NULL, *output_array = NULL;
	double *input_time   = NULL;
	double *input_pos[3] = { NULL }, *input_sig[3] = { NULL };

	float *x, *y;

	char testc[120];

	FILE *posfile = NULL;

	MATRIX *d = NULL;
	MATRIX *g = NULL;
	MATRIX *m = NULL;

/**/
	if ( argc < 2 ) {
		fprintf(stderr, "Usage: %s filename\n", argv[0]);
		return -1;
	}
/**/
	posfile = fopen(argv[1], "r");
	if ( posfile == NULL ) {
		fprintf(stderr, "Can't open the file %s\n", argv[1]);
		return -1;
	}
/**/
	while ( !feof(posfile) ) {
		fgets(testc, 119, posfile);
		count++;
	}
	rewind(posfile);

	input_array = (double *)malloc(sizeof(double) * count);
	input_time  = (double *)malloc(sizeof(double) * count);

	for ( i=0; i<3; i++ ) {
		input_pos[i] = (double *)malloc(sizeof(double) * count);
		input_sig[i] = (double *)malloc(sizeof(double) * count);
	}

	for ( i=0; i<count; i++ ) {
		fscanf(posfile, "%lf %lf %lf %lf %lf %lf %lf\n",
			input_time+i, input_pos[0]+i, input_sig[0]+i, input_pos[1]+i, input_sig[1]+i,
			input_pos[2]+i, input_sig[2]+i);

		if ( i == 0 ) {
			pivot[0] = *input_pos[0];
			pivot[1] = *input_pos[1];
			pivot[2] = *input_pos[2];
		}

		input_array[i]   = 1.0;
		input_pos[0][i] -= pivot[0];
		input_pos[1][i] -= pivot[1];
		input_pos[2][i] -= pivot[2];
	}

/**/
	x = (float *)malloc(sizeof(float)*count);
	y = (float *)malloc(sizeof(float)*count);
/* Just plot the result by PGPLOT library */
	if ( !cpgopen("position_with_time/cps") ) return -1;
	cpgsubp(1,3);

	for ( j=0; j<3; j++ ) {
	/* Assign the d matrix */
		d = matrix_new( count, 1 );
		matrix_assign_col( d, input_pos[j], 1, count );

	/* Assign the g matrix */
		g = matrix_new( count, NUM_PARAM );
	/* Parameter 0, Constant term */
		matrix_assign_col( g, input_array, 1, count );
	/* Parameter 1, Simple time term */
		matrix_assign_col( g, input_time, 2, count );
	/* Parameter 2, Seasoal shift sin term */
		matrix_assign_col( g, input_time, 3, count );
		matrix_apply_col( g, seasonal_term_sin, 3 );
	/* Parameter 3, Seasoal shift cos term */
		matrix_assign_col( g, input_time, 4, count );
		matrix_apply_col( g, seasonal_term_cos, 4 );
	/* Parameter 4, Co-seismic displacement term */
		matrix_assign_col( g, input_time, 5, count );
		matrix_apply_col( g, heavy_side_eq_jiashian, 5 );
	/* Parameter 5, Post-seismic displacement term */
		matrix_assign_col( g, input_time, 6, count );
		matrix_apply_col( g, post_seismic_def_eq_jiashian, 6 );

	/**/
		m = matrix_div( d, g );
		matrix_free( d );
		matrix_free( g );
		output_array = (double *)malloc(sizeof(double)*NUM_PARAM);
		matrix_extract_seq( m, output_array, NUM_PARAM );
		matrix_free( m );

	/**/
		for ( i=0; i<NUM_PARAM; i++)
			printf("%lf ", output_array[i]);

		printf("\n");

		max_y = min_y = input_pos[j][0];
		for ( i=0; i<count; i++ ) {
			x[i] = input_time[i];
			y[i] = input_pos[j][i];
			if ( y[i] > max_y )
				max_y = y[i];
			if ( y[i] < min_y )
				min_y = y[i];
		}
		pivot[0] = (max_y - min_y)*0.1;
		cpgsci(1);
		cpgenv( (float)input_time[0], (float)input_time[count - 1], min_y - (float)pivot[0], max_y + (float)pivot[0], 0, 0);
		cpglab("Year", "Position(mm)", "Position Change over Years in EW direc.");
		cpgpt(count, x, y, -1);
		for ( i=0; i<count; i++ ) {
			y[i] = final_function((double)x[i], output_array[0], output_array[1], output_array[2], output_array[3], output_array[4], output_array[5]);
		}
		cpgsci(2);
		cpgline(count, x, y);
		//sprintf(testc, "Y = %.3f + %.3f * X", output_array[0], output_array[1], output_array[2], output_array[3], output_array[4], output_array[5]);
		//cpgtext(2014.5, 10.0, testc);
		//sprintf(testc, "Velocity: %.3f cm/yr", output_array[1]);
		//cpgtext(2014.5, 9.5, testc);
	}
	cpgend();
	return 0;
}

static double seasonal_term_sin( const double input ) {
	return sin(2.0 * PI * input);
}

static double seasonal_term_cos( const double input ) {
	return cos(2.0 * PI * input);
}

static double heavy_side( const double input ) {
	return (input > 0.0 ? 1.0 : 0.0);
}

static double post_seismic_def( const double input ) {
	return exp(-input*TAU);
}

static double heavy_side_eq_jiashian( const double input ) {
	return heavy_side(input - JIASHIAN_EQ);
}

static double post_seismic_def_eq_jiashian( double input ) {
	input -= JIASHIAN_EQ;
	return heavy_side(input)*post_seismic_def(input);
}

static double final_function( const double x, const double a, const double b, const double c, const double d, const double e, const double f ) {
	return (a + b*x + c*seasonal_term_sin(x) + d*seasonal_term_cos(x) + e*heavy_side_eq_jiashian(x) + f*post_seismic_def_eq_jiashian(x));
}
