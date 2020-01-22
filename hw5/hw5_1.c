#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include <constant.h>
#include <formulas.h>
#include <matrix.h>
#include <okada91.h>
#include <faultgeom.h>

static FAULT_MODEL *faultmodel_gen( const FAULT_GEOM *, const double, const double, const double );

static int arg_parse( int, char ** );
static int stacoor_file_parse( const char *);
static int stadisp_file_parse( const char *);
static double misfit_compute( MATRIX *, MATRIX * );

static void end_process( void );

static double Length[3];
static double Depth[3];
static double Dip[3];
static double Strike[3];
static double East[3];
static double North[3];

static int       StaCount = 0;
static STA_COOR *StaCoor = NULL;
static double   *InputDisp = NULL;

/**/
static int arg_parse( int, char ** );

int main( int argc, char **argv )
{
	int     i;
	double  tmp, prev_misfit = 1e6;

	struct timespec TT, TT2;            /* Nanosecond Timer */

	MATRIX *d  = NULL;
	MATRIX *dhat = NULL;
	//MATRIX *wd = NULL;
	MATRIX *g  = NULL;
	MATRIX *s  = NULL;
	//MATRIX *w  = NULL;
	//MATRIX *wg = NULL;

	double length;
	double depth;
	double dip;
	double strike;
	double east;
	double north;


	FAULT_GEOM faultm = { 0.0 };
	FAULT_GEOM *faultp = NULL;

/**/
	if ( arg_parse( argc, argv ) ) return -1;

	fprintf(stdout, "Processing...\n");
	clock_gettime(CLOCK_REALTIME, &TT);

	double *us = (double *)malloc(sizeof(double)*StaCount*3);
	double *ud = (double *)malloc(sizeof(double)*StaCount*3);
	double dd[9];
	double ss[6];
/**/
/* Assign the weighted matrix */
/*
	w = matrix_new( DataCount, DataCount );
	matrix_assign_diag( w, InputSig[j], DataCount );
	matrix_apply_diag( w, self_invert );
*/
/* Assign the d matrix */
	d = matrix_new( StaCount*3, 1 );
	matrix_assign_seq( d, InputDisp, StaCount*3);

/* Assign the g matrix */
	g = matrix_new( StaCount*3, 2 );

	length = Length[0];

	for ( depth = Depth[0]; depth <= Depth[2]; depth += Depth[1] )
		for ( strike = Strike[0]; strike <= Strike[2]; strike += Strike[1] )
			for ( north = North[0]; north <= North[2]; north += North[1] )
				for ( east = East[0]; east <= East[2]; east += East[1] )
					for ( dip = Dip[0]; dip <= Dip[2]; dip += Dip[1] ) {
					/* Parameter 0, Constant term */
						fault2patch( transfault_bottom( (FAULT_GEOM){ length, depth/sin(dip*DEG2RAD), 0.0, -180.0+dip, strike, east, north } ), 1, 1, &faultp );
						//printf("%lf %lf %lf %lf %lf %lf %lf\n", faultp->length, faultp->width, faultp->depth, faultp->dip, faultp->strike, faultp->east, faultp->north);

						FAULT_MODEL *f1 = faultmodel_gen( faultp, 1.0, 0.0, 0.0 );
						FAULT_MODEL *f2 = faultmodel_gen( faultp, 0.0, 1.0, 0.0 );

						for ( i = 0; i < StaCount; i++ ) {
							disloc3d( f1 ,StaCoor+i, 1.0, 0.25, us+i*3, dd, ss );
							disloc3d( f2 ,StaCoor+i, 1.0, 0.25, ud+i*3, dd, ss );

							//printf("%.20lf %.20lf %.20lf\n", (us+i*3)[0], (us+i*3)[1], (us+i*3)[2]);
						}

						free(f1);
						free(f2);

						matrix_assign_col( g, us, 1, StaCount*3 );
						matrix_assign_col( g, ud, 2, StaCount*3 );

					/**/
						s = matrix_div( d, g );

						dhat = matrix_mul( g, s );
						tmp = misfit_compute( d, dhat );
						matrix_free( s );
						matrix_free( dhat );
						//fprintf(stdout, "Misfit: %lf\n", tmp);
						if ( tmp < prev_misfit ) {
							prev_misfit = tmp;
							faultm = (FAULT_GEOM){ length, depth/sin(dip*DEG2RAD), 0.0, -180.0+dip, strike, east, north };
						}
					}
/**/
	free(us);
	free(ud);
	free(faultp);
	matrix_free( d );
	matrix_free( g );
	end_process();

	fprintf(stdout, "Best result:\n");
	fprintf(stdout, "	%lf %lf %lf %lf %lf %lf %lf\n", faultm.length, faultm.width, faultm.depth, faultm.dip, faultm.strike, faultm.east, faultm.north);
	fprintf(stdout, "	Misfit: %lf\n", prev_misfit);


/* Nanosecond Timer */
	clock_gettime(CLOCK_REALTIME, &TT2);
	printf("Process Time: %.6lf sec\n", (TT2.tv_sec - TT.tv_sec) + (TT2.tv_nsec - TT.tv_nsec)*1.0e-9);

	return 0;
}

/**/
static int arg_parse( int argc, char **argv ) {
/**/

	if ( argc < 3 ) {
		fprintf(stderr, "Usage: %s station_coor_file station_disp_file\n", argv[0]);
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


static int stacoor_file_parse( const char *filename ) {
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
			StaCount++;
	}
	StaCount++;
	rewind(infile);

/**/
	StaCoor = (STA_COOR *)malloc(sizeof(STA_COOR) * StaCount);
/**/
	for ( i=0; i<StaCount; i++ ) {
		if ( fscanf(infile, "%lf %lf %lf\n", &((StaCoor + i)->x), &((StaCoor + i)->y), &((StaCoor + i)->z)) != 3 ) return -1;
	/**/
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

	return;
}
