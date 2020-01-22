#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include <constant.h>
#include <formulas.h>
#include <matrix.h>
#include <okada91.h>
#include <faultgeom.h>

#define  ORIGIN_LAT    23.0l
#define  ORIGIN_LON   120.7l

static FAULT_MODEL *faultmodel_gen( const FAULT_GEOM *, const double, const double, const double );

static double weighting_gen( const double );

static int arg_parse( int, char ** );
static int stacoor_file_parse( const char *);
static int stadisp_file_parse( const char *);
static int dispcov_file_parse( const char *);

static double misfit_compute( MATRIX *, MATRIX * );

static void end_process( void );

static int       StaCount = 0;
static STA_COOR *StaCoor = NULL;
static double   *InputDisp = NULL;
static double   *InputCov = NULL;
static int       Patch_d, Patch_s;
static double    Alpha = 1.0;
/**/
static int arg_parse( int, char ** );

int main( int argc, char **argv )
{
	int     i, j, npatches;
	double  tmp, prev_misfit = 1e6;

	struct timespec TT, TT2;            /* Nanosecond Timer */

	MATRIX *d  = NULL;
	MATRIX *dhat = NULL;
	MATRIX *wd = NULL;
	MATRIX *wwd = NULL;
	MATRIX *g  = NULL;
	MATRIX *s  = NULL;
	MATRIX *w  = NULL;
	MATRIX *wg = NULL;
	MATRIX *wwg = NULL;

	double length =  30.0;
	double tdepth =  12.0;
	double bdepth =  20.0;
	double dip    =  15.0;
	double strike = 275.0;
	double east   = -24.0;
	double north  = -12.0;
	double width  = (bdepth - tdepth)/sin(dip*DEG2RAD);

	FAULT_GEOM faultm = { 0.0 };
	FAULT_GEOM *faultp = NULL;

/**/
	if ( arg_parse( argc, argv ) ) return -1;

	fprintf(stdout, "Processing...\n");
	clock_gettime(CLOCK_REALTIME, &TT);

	double *laplacian = NULL, *tmparray = NULL;
	double *us = (double *)malloc(sizeof(double)*StaCount*3);
	double *ud = (double *)malloc(sizeof(double)*StaCount*3);
	double dd[9];
	double ss[6];
/**/
/* Assign the weighted matrix */
	w = matrix_new( StaCount*3, StaCount*3 );
	matrix_assign_diag( w, InputCov, StaCount*3 );
	matrix_apply_diag( w, weighting_gen );

/* Move fault coordinates ref. to the bottom, and divide it into subfaults */
	faultm = (FAULT_GEOM){ length, width, tdepth, -180.0+dip, strike, east, north };
	npatches = fault2patch( transfault_bottom( faultm ), Patch_s, Patch_d, &faultp );
	//printf("%lf %lf %lf %lf %lf %lf %lf\n", faultp->length, faultp->width, faultp->depth, faultp->dip, faultp->strike, faultp->east, faultp->north);

	laplacian = (double *)malloc(sizeof(double)*npatches*2);

/* Assign the d matrix */
	d = matrix_new( StaCount*3, 1 );
	matrix_assign_seq( d, InputDisp, StaCount*3);
	wd = matrix_mul( w, d );

/* Extend the matrix */
	wwd = matrix_new( StaCount*3 + npatches*2, 1 );
	tmparray = (double *)malloc(sizeof(double)*StaCount*3);
	matrix_extract_seq( wd, tmparray, StaCount*3 );
	matrix_assign_seq( wwd, tmparray, StaCount*3 );
	free(tmparray);
	matrix_free( wd );

	g = matrix_new( StaCount*3, npatches*2 );
/* Generating Green's function */
	for ( j = 0; j < npatches; j++ ) {
	/**/
		FAULT_MODEL *f1 = faultmodel_gen( faultp, 1.0, 0.0, 0.0 );
		FAULT_MODEL *f2 = faultmodel_gen( faultp, 0.0, 1.0, 0.0 );
	/**/
		for ( i = 0; i < StaCount; i++ ) {
			disloc3d( f1 ,StaCoor+i, 1.0, 0.25, us+i*3, dd, ss );
			disloc3d( f2 ,StaCoor+i, 1.0, 0.25, ud+i*3, dd, ss );
			//printf("%le %le %.le\n", (us+i*3)[0], (us+i*3)[1], (us+i*3)[2]);
		}
	/**/
		free(f1);
		free(f2);
	/* Assign the g matrix */
		matrix_assign_col( g, us, j+1, StaCount*3 );
		matrix_assign_col( g, ud, npatches+j+1, StaCount*3 );
	}

/* Apply weighting matrix */
	wg = matrix_mul( w, g );

/* Extend the matrix */
	wwg = matrix_new( StaCount*3 + npatches*2, npatches*2  );
	tmparray = (double *)malloc(sizeof(double)*StaCount*npatches*6);
	matrix_extract_seq( wg, tmparray, StaCount*npatches*6 );
	matrix_assign_seq( wwg, tmparray, StaCount*npatches*6 );
	free(tmparray);
	matrix_free( wg );

	npatches = npatches > 1 ? npatches : 0;
/* Generating the smoothing Laplacian matrix */
	for ( j = 0; j < npatches; j++ ) {
		laplacian[j] = -2.0 * Alpha;
	/* Along Strike */
		for ( i = 0; i < npatches; i++ ) {
			if ( i != j ) {
				if ( i == (j - Patch_d) || i == (j + Patch_d) )
					laplacian[i] = 1.0 * Alpha;
				else
					laplacian[i] = 0.0;
			}
		}
		matrix_assign_row( wwg, laplacian, StaCount*3+j+1, npatches );
		memset(laplacian, 0, npatches*sizeof(double));

	/* Along Dip */
		laplacian[npatches+j] = -2.0 * Alpha;
		for ( i = 0; i < npatches; i++ ) {
			if ( i != j ) {
				if ( (i % Patch_d) == (j % Patch_d) && (i == (j - 1) || i == (j + 1)) )
					laplacian[npatches+i] = 1.0 * Alpha;
				else
					laplacian[npatches+i] = 0.0;
			}
		}
		matrix_assign_row( wwg, laplacian, StaCount*3+npatches+j+1, npatches*2 );
	}
/*
	Inversion for fault slip, and the left-lateral slip & thrust slip
	is positive so we should use Non-negative Least-squares here.
*/
	//s = matrix_div( d, g );
	s = matrix_nnls( wwd, wwg );

	dhat = matrix_mul( g, s );
	tmp = misfit_compute( d, dhat );
	matrix_free( s );
	matrix_free( dhat );
	fprintf(stdout, "Misfit: %lf\n", tmp);

/**/
	free(us);
	free(ud);
	free(faultp);
	matrix_free( d );
	matrix_free( g );
	matrix_free( wwd );
	matrix_free( wwg );
	end_process();

/* Nanosecond Timer */
	clock_gettime(CLOCK_REALTIME, &TT2);
	printf("Process Time: %.6lf sec\n", (TT2.tv_sec - TT.tv_sec) + (TT2.tv_nsec - TT.tv_nsec)*1.0e-9);

	return 0;
}

/**/
static int arg_parse( int argc, char **argv ) {
/**/

	if ( argc != 7 ) {
		fprintf(stderr, "Usage: %s station_coor_file station_disp_file dispcov_file patch_num_dip patch_num_strike alpha\n", argv[0]);
		return -1;
	}

/**/
	if ( stacoor_file_parse( argv[1] ) ) {
		fprintf(stderr, "File parsing error!!\n");
		return -1;
	}

/**/
	if ( stadisp_file_parse( argv[2] ) ) {
		fprintf(stderr, "File parsing error!!\n");
		return -1;
	}

/**/
	if ( dispcov_file_parse( argv[3] ) ) {
		fprintf(stderr, "File parsing error!!\n");
		return -1;
	}

	Patch_d = atoi(argv[4]);
	Patch_s = atoi(argv[5]);
	Alpha   = atof(argv[6]);

/*
	fprintf(stdout, "Please input the searching fault length:\n");
	if ( fscanf(stdin, "%lf", &Length[0]) != 1 ) return -1;

	fprintf(stdout, "Please input the searching range of fault bottom depth(lbound inc ubound):\n");
	if ( fscanf(stdin, "%lf %lf %lf", &Depth[0], &Depth[1], &Depth[2]) < 3 ) return -1;

	fprintf(stdout, "Please input the searching range of fault dip angle(lbound inc ubound):\n");
	if ( fscanf(stdin, "%lf %lf %lf", &Dip[0], &Dip[1], &Dip[2]) < 3 ) return -1;

	fprintf(stdout, "Please input the searching range of fault strike angle(lbound inc ubound):\n");
	if ( fscanf(stdin, "%lf %lf %lf", &Strike[0], &Strike[1], &Strike[2]) < 3 ) return -1;

	fprintf(stdout, "Please input the searching range of fault x coordinate of midpoint(lbound inc ubound):\n");
	if ( fscanf(stdin, "%lf %lf %lf", &East[0], &East[1], &East[2]) < 3 ) return -1;

	fprintf(stdout, "Please input the searching range of fault y coordinates of midpoint(lbound inc ubound):\n");
	if ( fscanf(stdin, "%lf %lf %lf", &North[0], &North[1], &North[2]) < 3 ) return -1;
*/


	return 0;
}

static FAULT_MODEL *faultmodel_gen( const FAULT_GEOM *faultg, const double disl1, const double disl2, const double disl3 )
{
	FAULT_MODEL *ret = (FAULT_MODEL *)malloc(sizeof(FAULT_MODEL));

	ret->length = faultg->length;
	ret->width  = faultg->width;
	ret->depth  = faultg->depth;
	ret->dip    = faultg->dip;
	ret->strike = faultg->strike;
	ret->east   = faultg->east;
	ret->north  = faultg->north;
	ret->disl1  = disl1;
	ret->disl2  = disl2;
	ret->disl3  = disl3;

	return ret;
}

/*
*/
static double weighting_gen( const double input ) {
	if ( fabs(input) > DBL_EPSILON )
		return 1.0 / sqrt(input);
	else
		return 0.0;
}

/*
*/
static int stacoor_file_parse( const char *filename ) {
	int   i;
	char  c;
	double _lat, _lon, _z;
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
			StaCount++;
	}
	StaCount++;
	rewind(infile);

/**/
	StaCoor = (STA_COOR *)malloc(sizeof(STA_COOR) * StaCount);
/**/
	for ( i=0; i<StaCount; i++ ) {
		if ( fscanf(infile, "%lf %lf %lf\n", &_lon, &_lat, &_z) != 3 ) return -1;

		(StaCoor + i)->x = coor2distf( ORIGIN_LAT, ORIGIN_LON, ORIGIN_LAT, _lon ) * 1e-5 * ( _lon > ORIGIN_LON ? 1.0 : -1.0 ); /* unit: km */
		(StaCoor + i)->y = coor2distf( ORIGIN_LAT, ORIGIN_LON, _lat, ORIGIN_LON ) * 1e-5 * ( _lat > ORIGIN_LAT ? 1.0 : -1.0 ); /* unit: km */
	/**/(StaCoor + i)->z = _z;
	}
/**/
	fclose(infile);
	return 0;
}

static int stadisp_file_parse( const char *filename ) {
	int    i;
	FILE   *infile;
	double *inputptr = NULL;

/**/
	infile = fopen(filename, "r");
	if ( infile == NULL ) {
		fprintf(stderr, "Can not open the file %s!\n", filename);
		return -1;
	}

/*
	for ( c=getc(infile); c!=EOF; c=getc(infile) ) {
		if ( c == '\n' )
			StaCount++;
	}
	rewind(infile);
*/
/**/
	InputDisp = (double *)malloc(sizeof(double) * StaCount * 3);
/**/
	inputptr = InputDisp;
	for ( i = 0; i < StaCount; i++ ) {
		if ( fscanf(infile, "%lf %lf %lf\n", inputptr + i*3, inputptr + i*3 + 1, inputptr + i*3 + 2) != 3 ) return -1;
	/**/
	}
/**/
	fclose(infile);
	return 0;
}

/*
*/
static int dispcov_file_parse( const char *filename ) {
	int    i;
	FILE   *infile;
	double *inputptr = NULL;

/**/
	infile = fopen(filename, "r");
	if ( infile == NULL ) {
		fprintf(stderr, "Can not open the file %s!\n", filename);
		return -1;
	}

/*
	for ( c=getc(infile); c!=EOF; c=getc(infile) ) {
		if ( c == '\n' )
			StaCount++;
	}
	rewind(infile);
*/
/**/
	InputCov = (double *)malloc(sizeof(double) * StaCount * 3);
/**/
	inputptr = InputCov;
	for ( i = 0; i < StaCount; i++ ) {
		if ( fscanf(infile, "%lf %lf %lf\n", inputptr + i*3, inputptr + i*3 + 1, inputptr + i*3 + 2) != 3 ) return -1;
	/**/
	}
/**/
	fclose(infile);
	return 0;
}

static double misfit_compute( MATRIX *obs, MATRIX* pred ) {
	int     i;
	double  res      = 0.0;
	double *misfit_a = (double *)malloc(sizeof(double)*StaCount*3);
	MATRIX *misfit_m = matrix_sub( obs, pred );

	//matrix_apply_col( misfit_m, self_square, 1 );
	matrix_extract_seq( misfit_m, misfit_a, StaCount*3 );
	for ( i=0; i<StaCount*3; i++ ) res += fabs(misfit_a[i]);
	res /= (double)(StaCount*3.0);
	res *= 1000.0;

	free(misfit_a);
	matrix_free( misfit_m );

	return res;
}

static void end_process( void ) {
	free(StaCoor);
	free(InputDisp);
	free(InputCov);

	return;
}
