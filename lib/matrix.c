#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <matrix.h>


static int matrix_same_size( const MATRIX *, const MATRIX * );
static int matrix_square( const MATRIX * );
static int matrix_copy( MATRIX *, const MATRIX * );

/**/
MATRIX *matrix_new( const int row, const int column ) {
	MATRIX *res;

	res          = (MATRIX *)calloc(1, sizeof(MATRIX));
	res->i       = row;
	res->j       = column;
	res->total   = row * column;
	res->element = (double *)calloc(res->total, sizeof(double));

	return res;
}

/**/
MATRIX *matrix_identity( const int rank ) {
	MATRIX *res = matrix_new( rank, rank );
	int     i;

	for ( i=0; i<rank; i++ ) res->element[i * rank + i] = 1.0;

	return res;
}

/**/
MATRIX *matrix_add( const MATRIX *a, const MATRIX *b ) {
	int     i, j;
	MATRIX *res = NULL;

	if ( matrix_same_size(a, b) ) {
		res = matrix_new( a->i, a->j );

		for ( i=0; i<a->i; i++ ) {
			for ( j=0; j<a->j; j++ )
				res->element[i * a->j + j] = a->element[i * a->j + j] + b->element[i * a->j + j];
		}
	}

	return res;
}

/**/
MATRIX *matrix_sub( const MATRIX *a, const MATRIX *b ) {
	int     i, j;
	MATRIX *res = NULL;

	if ( matrix_same_size(a, b) ) {
		res = matrix_new( a->i, a->j );

		for ( i=0; i<a->i; i++ ) {
			for ( j=0; j<a->j; j++ )
				res->element[i * a->j + j] = a->element[i * a->j + j] - b->element[i * a->j + j];
		}
	}

	return res;
}

/**/
MATRIX *matrix_mul( const MATRIX *a, const MATRIX *b ) {
	int     i, j, k;
	MATRIX *res = NULL;

	if ( a->j == b->i ) {
		res = matrix_new( a->i, b->j );

		for ( i=0; i<a->i; i++ ) {
			for ( j=0; j<b->j; j++ ) {
				res->element[i*(b->j) + j] = 0.0;
				for ( k=0; k<a->j; k++ ) {
					res->element[i * b->j + j] += a->element[i * a->j + k] * b->element[k * b->j + j];
				}
			}
		}
	}

	return res;
}

/**/
MATRIX *matrix_transpose( const MATRIX *a ) {
	int     i, j;
	MATRIX *res = NULL;

	res = matrix_new( a->j, a->i );

	for ( i=0; i<a->i; i++ ) {
		for ( j=0; j<a->j; j++ ) {
			res->element[j * a->i + i] = a->element[i * a->j + j];
		}
	}

	return res;
}

/*
MATRIX *matrix_adjugate( const MATRIX *a ) {
	int     i, j;
	MATRIX *res = NULL;

	if ( matrix_square( a ) ) {
		res = matrix_init( a->i, a->j );
		for ( i=0; i<a->i; i++ ) {
			for ( j=0; j<a->j; j++ ) {
				res->element[j*(a->i) + i] = a->element[i*(a->j) + j];
			}
		}
	}

	return res;
}
*/
/**/
MATRIX *matrix_inverse( const MATRIX *a ) {
	int     i, j, k;
	int     prow = 0;
	double  pivot;
	MATRIX *res = NULL;
	MATRIX *tmp = NULL;

	if ( matrix_square( a ) ) {
		tmp = matrix_new( a->i, a->j );
		res = matrix_identity( a->i );
		matrix_copy( tmp, a );

		for ( i=0; i<a->i; i++ ) {
		/* */
			pivot = tmp->element[i * a->j + i];
			for ( j=i, prow=i; j<a->i; j++ ) {
				if ( tmp->element[j * a->j + i] > pivot ) {
					pivot = tmp->element[j * a->j + i];
					prow  = j;
				}
			}
		/* */
			if ( prow != i ) {
				for ( j=0; j<a->j; j++ ) {
				/* */
					pivot                         = tmp->element[i * a->j + j];
					tmp->element[i * a->j + j]    = tmp->element[prow * a->j + j];
					tmp->element[prow * a->j + j] = pivot;
				/* */
					pivot                         = res->element[i * a->j + j];
					res->element[i * a->j + j]    = res->element[prow * a->j + j];
					res->element[prow * a->j + j] = pivot;
				}
			}
		/* */
			pivot = tmp->element[i * a->j + i];
			for ( j=0; j<a->j; j++ ) {
			/* */
				tmp->element[i * a->j + j] /= pivot;
			/* */
				res->element[i * a->j + j] /= pivot;
			}
		/* */
			for ( j=i+1; j<a->i; j++ ) {
				pivot = tmp->element[j * a->j + i];
				if ( fabs(pivot) > DBL_EPSILON ) {
					for ( k=0; k<a->j; k++ ) {
					/* */
						if ( fabs(tmp->element[i * a->j + k]) > DBL_EPSILON )
							tmp->element[j * a->j + k] -= pivot * tmp->element[i * a->j + k];
					/* */
						if ( fabs(res->element[i * a->j + k]) > DBL_EPSILON )
							res->element[j * a->j + k] -= pivot * res->element[i * a->j + k];
					}
				}
			}
		}

		prow = a->i - 1;
		for ( i=0; i<prow; i++ ) {
			for ( j=i+1; j<a->i; j++ ) {
				pivot = tmp->element[i * a->j + j];
				if ( fabs(pivot) > DBL_EPSILON ) {
					for ( k=0; k<a->j; k++ ) {
					/* */
						if ( fabs(tmp->element[j * a->j + k]) > DBL_EPSILON )
							tmp->element[i * a->j + k] -= pivot * tmp->element[j * a->j + k];
					/* */
						if ( fabs(res->element[j * a->j + k]) > DBL_EPSILON )
							res->element[i * a->j + k] -= pivot * res->element[j* a->j + k];
					}
				}
			}
		}

		matrix_free( tmp );
	}

	return res;
}

MATRIX *matrix_div( const MATRIX *a, const MATRIX *b ) {
	MATRIX *gt     = matrix_transpose( b );
	MATRIX *gtg    = matrix_mul( gt, b );
	MATRIX *igtg   = matrix_inverse( gtg );
	MATRIX *igtggt = matrix_mul( igtg, gt );
	MATRIX *res    = matrix_mul( igtggt, a );

	matrix_free( gt );
	matrix_free( gtg );
	matrix_free( igtg );
	matrix_free( igtggt );

	return res;
}

/* Assignment functions */
MATRIX *matrix_assign_seq( MATRIX *dest, const double *src, const int data_size ) {
	if ( dest->total > data_size ) {
		memcpy(dest->element, src, data_size * sizeof(double));
		return dest;
	}

	return NULL;
}

MATRIX *matrix_assign_row( MATRIX *dest, const double *src, int row_index, const int data_size ) {
	int j;

	row_index--;

	if ( dest->j == data_size ) {
		memcpy(dest->element + row_index * dest->j, src, dest->j * sizeof(double));
		return dest;
	}
	else if ( dest->j > data_size ) {
		memcpy(dest->element + row_index * dest->j, src, data_size * sizeof(double));
		for ( j=data_size; j<dest->j; j++ )
			dest->element[row_index * dest->j + j] = 0.0;
		return dest;
	}

	return NULL;
}

MATRIX *matrix_assign_col( MATRIX *dest, const double *src, int col_index, const int data_size ) {
	int i;

	col_index--;

	if ( dest->i == data_size ) {
		for ( i=0; i<dest->i; i++ )
			dest->element[i * dest->j + col_index] = src[i];
		return dest;
	}
	else if ( dest->i > data_size ) {
		for ( i=0; i<dest->i; i++ ) {
			if ( i < data_size )
				dest->element[i * dest->j + col_index] = src[i];
			else
				dest->element[i * dest->j + col_index] = 0.0;
		}
		return dest;
	}

	return NULL;
}

MATRIX *matrix_assign_diag( MATRIX *dest, const double *src, const int data_size ) {
	int i;

	if ( matrix_square( dest ) ) {
		if ( dest->i == data_size ) {
			for ( i=0; i<dest->i; i++ )
				dest->element[i * dest->j + i] = src[i];
			return dest;
		}
		else if ( dest->i > data_size ) {
			for ( i=0; i<dest->i; i++ ) {
				if ( i < data_size )
					dest->element[i * dest->j + i] = src[i];
				else
					dest->element[i * dest->j + i] = 0.0;
			}
			return dest;
		}
	}

	return NULL;
}

double *matrix_extract_seq( const MATRIX *src, double *dest, const int dest_size ) {
	if ( src->total <= dest_size ) {
		memcpy(dest, src->element, dest_size * sizeof(double));
		return dest;
	}

	return NULL;
}

double matrix_determinant( const MATRIX *a ) {
	if ( matrix_square( a ) ) {
		return 0;
	}
	return -1;
}

void matrix_free( MATRIX *a ) {
	free(a->element);
	free(a);
	return;
}

static int matrix_same_size( const MATRIX *a, const MATRIX *b ) {
	return ( a->i == b->i && a->j == b->j );
}

static int matrix_square( const MATRIX *a ) {
	return ( a->i == a->j );
}

static int matrix_copy( MATRIX *dest, const MATRIX *src ) {
	if ( matrix_same_size( dest, src ) ) {
		memcpy(dest->element, src->element, src->total * sizeof(double));
		return 0;
	}
	return -1;
}