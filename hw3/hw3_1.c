#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>

/* PGPLOT C library header */
#include <cpgplot.h>
#include <constant.h>
#include <formulas.h>
#include <matrix.h>

#define DEFAULT_TAU  1.5f
#define BASIC_PARAM  4
#define JIASHIAN_EQ  2010.1726f
#define MEINONG_EQ   2016.1014f

/**/
static int arg_parse( int, char ** );
static int file_parse( const char *);

static float *array_d2f( float *, const double *, const int );
static double heavy_side( const double );
static double def_decay( const double );
static double self_invert( const double );
static double self_square( const double );

static double seasonal_term_sin( const double );
static double seasonal_term_cos( const double );
static double coseismic_def( const double );
static double postseismic_def( const double );
static double misfit_compute( MATRIX *, MATRIX * );

static void end_process( void );

static volatile double  Tau       = DEFAULT_TAU;
static volatile double *EqTimePtr = NULL;
static volatile int     EqCount   = 0;

static int     DataCount   = 0;
static double *InputEqTime = NULL;
static double *InputConst  = NULL;
static double *InputTime   = NULL;
static double *InputPos[3] = { NULL };
static double *InputSig[3] = { NULL };

struct timespec TT, TT2;            /* Nanosecond Timer */

int main( int argc, char **argv )
{
	int     i, j, iter;
	int     param_index = 0;
	int     total_param = BASIC_PARAM;
	float   min_y, max_y;
	float  *plot_x, *plot_y;
	double  tmp, prev_misfit = 0.0;
 	double *output_array = NULL;
	char    tmpstring[256];
	char   *dire_str[] = { "EW", "NS", "UD" };

	MATRIX *d  = NULL;
	MATRIX *wd = NULL;
	MATRIX *dp = NULL;
	MATRIX *g  = NULL;
	MATRIX *m  = NULL;
	MATRIX *w  = NULL;
	MATRIX *wg = NULL;

/**/
	if ( arg_parse( argc, argv ) ) return -1;
	total_param += 2 * EqCount;
	clock_gettime(CLOCK_REALTIME, &TT);
	//printf("There are total %d parameters!\n", total_param);
/**/
	output_array = (double *)malloc(sizeof(double)*DataCount);
	plot_x = (float *)malloc(sizeof(float)*DataCount);
	plot_y = (float *)malloc(sizeof(float)*DataCount);
/* Just plot the result by PGPLOT library */
	if ( !cpgopen("position_with_time/cps") ) return -1;
	cpgsubp(1,3);

	for ( j=0; j<3; j++ ) {
	/* Assign the weighted matrix */
		w = matrix_new( DataCount, DataCount );
		matrix_assign_diag( w, InputSig[j], DataCount );
		matrix_apply_diag( w, self_invert );

	/* Assign the d matrix */
		d = matrix_new( DataCount, 1 );
		matrix_assign_col( d, InputPos[j], 1, DataCount );
		wd = matrix_mul( w, d );

	/* Assign the g matrix */
		g = matrix_new( DataCount, total_param );
	/* Parameter 0, Constant term */
		matrix_assign_col( g, InputConst, 1, DataCount );
	/* Parameter 1, Simple time term */
		matrix_assign_col( g, InputTime, 2, DataCount );
	/* Parameter 2, Seasoal shift sin term */
		matrix_assign_col( g, InputTime, 3, DataCount );
		matrix_apply_col( g, seasonal_term_sin, 3 );
	/* Parameter 3, Seasoal shift cos term */
		matrix_assign_col( g, InputTime, 4, DataCount );
		matrix_apply_col( g, seasonal_term_cos, 4 );

		param_index = 5;
		for ( i=0; i<EqCount; i++ ) {
			EqTimePtr = InputEqTime + i;
		/* Parameter for Co-seismic displacement term */
			matrix_assign_col( g, InputTime, param_index, DataCount );
			matrix_apply_col( g, coseismic_def, param_index );
			param_index++;
		/* Parameter for Post-seismic displacement term */
			matrix_assign_col( g, InputTime, param_index, DataCount );
			matrix_apply_col( g, postseismic_def, param_index );
			param_index++;
		}

		iter = 0;
		Tau = DEFAULT_TAU;
		do {
			wg = matrix_mul( w, g );
		/**/
			m = matrix_div( wd, wg );
			matrix_free( wg );

			dp = matrix_mul( g, m );
			tmp = misfit_compute( d, dp );
			if ( fabs(tmp - prev_misfit) > (prev_misfit * 0.000001) && iter++ < 100 ) {
				matrix_free( dp );
				matrix_free( m );
				prev_misfit = tmp;
				Tau += 0.01;
				param_index = 6;
				for ( i=0; i<EqCount; i++ ) {
					EqTimePtr = InputEqTime + i;
				/* Parameter for Post-seismic displacement term */
					matrix_assign_col( g, InputTime, param_index, DataCount );
					matrix_apply_col( g, postseismic_def, param_index );
					param_index += 2;
				}
			}
			else break;
		} while ( 1 );
	/**/
		//printf("Tau = %lf\n", Tau);
		matrix_free( wd );
		matrix_free( w );
		matrix_free( d );
	/**/
		max_y = min_y = InputPos[j][0];
		for ( i=0; i<DataCount; i++ ) {
			plot_x[i] = InputTime[i];
			plot_y[i] = InputPos[j][i];
			if ( plot_y[i] > max_y ) max_y = plot_y[i];
			if ( plot_y[i] < min_y ) min_y = plot_y[i];
		}
		tmp = (max_y - min_y)*0.1;
	/**/
		cpgsci(1);
		cpgenv((float)InputTime[0], (float)InputTime[DataCount-1], min_y - (float)tmp, max_y + (float)tmp, 0, 0);
		sprintf(tmpstring, "Position Change over Years From %s", argv[1]);
		cpglab( j == 2 ? "Year" : "", "Position(mm)", j == 0 ? tmpstring : "" );
		cpgslw(2);
		cpgpt(DataCount, plot_x, plot_y, -1);
		cpgslw(1);
	/**/
		matrix_extract_seq( dp, output_array, DataCount );
		array_d2f( plot_y, output_array, DataCount );
		matrix_free( dp );
		matrix_free( g );
		cpgsci(2);
		cpgline(DataCount, plot_x, plot_y);

	/**/
		matrix_extract_seq( m, output_array, total_param );
		matrix_free( m );

		cpgsci(1);
		sprintf(tmpstring, "Direction %s", dire_str[j]);
		prev_misfit = tmp;
		tmp = (float)InputTime[0] + 0.1;
		cpgtext((float)tmp, (float)min_y + fabs(prev_misfit), tmpstring);
		sprintf(tmpstring, "Y = %.3f + (%.3f)X + (%.3f)sin(2*pi*X) + (%.3f)cos(2*pi*X) + (%.3f)H(X-Teq) + (%.3f)(1-exp((X-Teq)/(%.2f)))H(X-Teq)",
			output_array[0], output_array[1], output_array[2], output_array[3], output_array[4], output_array[5], Tau);
		cpgtext((float)tmp, (float)min_y + fabs(prev_misfit * 0.5), tmpstring);
		sprintf(tmpstring, "Velocity: %.3f mm/yr", output_array[1]);
		cpgtext(tmp, (float)min_y, tmpstring);
	/**/
		//for ( i=0; i<total_param; i++)
			//printf("%lf ", output_array[i]);
		//printf("\n");
	}
	cpgend();

	free(output_array);
	free(plot_x);
	free(plot_y);
	end_process();

/* Nanosecond Timer */
	clock_gettime(CLOCK_REALTIME, &TT2);
	printf("Process Time: %.6lf sec\n", (TT2.tv_sec - TT.tv_sec) + (TT2.tv_nsec - TT.tv_nsec)*1.0e-9);

	return 0;
}

/*
*/
static float *array_d2f( float *dest, const double *src, const int size ) {
	int i;
	for ( i=0; i<size; i++ ) dest[i] = (float)src[i];
	return dest;
}

/*
*/
static double heavy_side( const double input ) {
	return (input > 0.0 ? 1.0 : 0.0);
}

/*
*/
static double def_decay( const double input ) {
	return (1.0 - exp(-input*Tau));
}

/*
*/
static double self_invert( const double input ) {
	if ( fabs(input) > DBL_EPSILON )
		return 1.0 / input;
	else
		return 0.0;
}

/*
*/
static double self_square( const double input ) {
	return input * input;
}


/*
*/
static double seasonal_term_sin( const double input ) {
	return sin(2.0 * PI * input);
}

/*
*/
static double seasonal_term_cos( const double input ) {
	return cos(2.0 * PI * input);
}

/*
*/
static double coseismic_def( const double input ) {
	return heavy_side(input - *EqTimePtr);
}

/*
*/
static double postseismic_def( const double input ) {
	double _input = input - *EqTimePtr;
	return heavy_side(_input)*def_decay(_input);
}

/**/
static int arg_parse( int argc, char **argv ) {
	int    i;
	int    _eqDataCount = 0;
	double _eqtime  = 0.0;

/**/
	if ( argc < 2 ) {
		fprintf(stderr, "Usage: %s inputfile [#eq eq1_time eq2_time...]\n", argv[0]);
		return -1;
	}

/**/
	if ( file_parse( argv[1] ) ) {
		fprintf(stderr, "File parsing error!!\n");
		return -1;
	}

/**/
	if ( argc > 2 ) {
		EqCount = atoi(argv[2]);
		if ( argc != (3 + EqCount) ) {
			fprintf(stderr, "Number of Eq. doesn't match the following Eq. times!\n");
			return -1;
		}

		InputEqTime = (double *)malloc(sizeof(double)*EqCount);
		for ( i=0; i<EqCount; i++ ) {
			_eqtime = atof(argv[3+i]);
			if ( _eqtime > InputTime[0] && _eqtime < InputTime[DataCount-1] ) {
				InputEqTime[_eqDataCount] = _eqtime;
				_eqDataCount++;
			}
		}
		EqCount = _eqDataCount;
	}

	return 0;
}

static int file_parse( const char *filename ) {
	int   i;
	char  c;
	FILE *infile;

/**/
	infile = fopen(filename, "r");
	if ( infile == NULL ) {
		fprintf(stderr, "Can not open the file %s!\n", filename);
		return -1;
	}

/**/
	for ( c=getc(infile); c!=EOF; c=getc(infile) ) {
		if ( c == '\n' )
			DataCount++;
	}
	rewind(infile);

/**/
	InputConst = (double *)malloc(sizeof(double) * DataCount);
	InputTime  = (double *)malloc(sizeof(double) * DataCount);
	for ( i=0; i<3; i++ ) {
		InputPos[i] = (double *)malloc(sizeof(double) * DataCount);
		InputSig[i] = (double *)malloc(sizeof(double) * DataCount);
	}
/**/
	for ( i=0; i<DataCount; i++ ) {
		if ( fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf\n", InputTime+i, InputPos[0]+i, InputSig[0]+i,
				InputPos[1]+i, InputSig[1]+i, InputPos[2]+i, InputSig[2]+i) != 7 ) return -1;
	/**/
		InputConst[i] = 1.0;
		if ( i > 0 ) {
			InputPos[0][i] -= InputPos[0][0];
			InputPos[1][i] -= InputPos[1][0];
			InputPos[2][i] -= InputPos[2][0];
		}
	}
/**/
	InputPos[0][0] = 0.0;
	InputPos[1][0] = 0.0;
	InputPos[2][0] = 0.0;

	fclose(infile);
	return 0;
}

static double misfit_compute( MATRIX *obs, MATRIX* pred ) {
	int     i;
	double  res      = 0.0;
	double *misfit_a = (double *)malloc(sizeof(double)*DataCount);
	MATRIX *misfit_m = matrix_sub( obs, pred );

	matrix_apply_col( misfit_m, self_square, 1 );
	matrix_extract_seq( misfit_m, misfit_a, DataCount );
	for ( i=0; i<DataCount; i++ ) res += misfit_a[i];

	free(misfit_a);
	matrix_free( misfit_m );

	return res;
}

static void end_process( void ) {
	int i;

	free(InputEqTime);
	free(InputConst);
	free(InputTime);

	for ( i=0; i<3; i++ ) {
		free(InputPos[i]);
		free(InputSig[i]);
	}

	return;
}
