#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "../include/constant.h"
#include "../include/okada91.h"

#define L_EPSILON  1e-6

/**/
static void dccon0(double, double);
static void dccon2(double *, double *, double *, double, double);
static void ua(double, double, double, double, double, double, double *);
static void ub(double, double, double, double, double, double, double *);
static void uc(double, double, double, double, double, double, double, double *);

/* C0 */
static double ALP[5];
static double SD, CD;
static double SDSD, CDCD;
static double SDCD, S2D, C2D;

/* C1 */
static double XI2, ET2, Q2;
static double R, R2, R3, R5;
static double Y, D, TT;
static double ALX, ALE;
static double X11, Y11, X32, Y32;
static double EY, EZ, FY, FZ, GY, GZ, HY, HZ;

/**/
int disloc3d (
	FAULT_MODEL model,
	STA_COOR station,
	double mu,
	double nu,
	double u[3],
	double d[9],
	double s[6]
) {
	const double lambda  = mu * nu / (0.5 - nu); /* Origin: 2.0*mu*nu/(1.0-2.0*nu) */
	const double alpha   = (lambda + mu) / (lambda + 2.0 * mu);
	const double strike  = (model.strike - 90.0) * DEG2RAD;
	const double cos_s   = cos(strike);
	const double sin_s   = sin(strike);
	const double cos_d   = cos(model.dip * DEG2RAD);
	const double sin_d   = sin(model.dip * DEG2RAD);
	const double al      = model.length * 0.5;
	const double aw      = model.width * 0.5;
	const double depth   = model.depth - aw * sin_d;
	const double x       = cos_s * (-model.east + station.x) - sin_s * (-model.north + station.y);
	const double y       = -0.5 * cos_d * model.width + sin_s * (-model.east + station.x) + cos_s * (-model.north + station.y);

	double uu[12] = { 0.0 };

	if ( (model.depth - sin_d * model.width) < -L_EPSILON ||
		model.length < DBL_EPSILON ||
		model.width < DBL_EPSILON ||
		model.depth < 0.0 )
	{
		fprintf(stderr, "disloc3d: Unphysical model!\n");
		return -1;
	}

/* */
	if ( dc3d( alpha, x, y, station.z, depth, model.dip, al, al, aw, aw, model.disl1, model.disl2, model.disl3,
		&uu[0], &uu[1], &uu[2], &uu[3], &uu[4], &uu[5], &uu[6], &uu[7], &uu[8], &uu[9], &uu[10], &uu[11] ) );
		//fprintf(stderr, "disloc3d: Singular result!\n");

/* */
	u[0] = cos_s * uu[0] + sin_s * uu[1];
	u[1] = -sin_s * uu[0] + cos_s * uu[1];
	u[2] = uu[2];

	//printf("%le %le %le\n", u[0], u[1], u[2]);

/**/
	d[0] = cos_s * cos_s * uu[3] + cos_s * sin_s * (uu[4] + uu[6]) + sin_s * sin_s * uu[7];
	d[1] = cos_s * cos_s * uu[4] - sin_s * sin_s * uu[6] + cos_s * sin_s * (-uu[3] + uu[7]);
	d[2] = cos_s * uu[5] + sin_s * uu[8];
	d[3] = -(sin_s * (cos_s * uu[3] + sin_s * uu[4])) + cos_s * (cos_s * uu[6] + sin_s * uu[7]);
	d[4] = sin_s * sin_s * uu[3] - cos_s * sin_s * (uu[4] + uu[6]) + cos_s * cos_s * uu[7];
	d[5] = -(sin_s * uu[5]) + cos_s * uu[8];
	d[6] = cos_s * uu[9] + sin_s * uu[10];
	d[7] = -(sin_s * uu[9]) + cos_s * uu[10];
	d[8] = uu[11];

/**/
	double theta = lambda * (d[0] + d[4] + d[8]);

	s[0] = theta + 2.0 * mu * d[0];
	s[1] = mu * (d[1] + d[3]);
	s[2] = mu * (d[2] + d[6]);
	s[3] = theta + 2.0 * mu * d[4];
	s[4] = mu * (d[5] + d[7]);
	s[5] = theta + 2.0 * mu * d[8];

	return 0;
}


/**/
int dc3d (
	const double _alpha,
	double _x, double _y, double z,
	double _depth, double _dip, double _al1, double _al2, double _aw1, double _aw2,
	double _disl1, double _disl2, double _disl3,
	double u[12]
/*
	*_ux, *_uy, *_uz, *_uxx, *_uyx, *_uzx,
	*_uxy, *_uyy, *_uzy, *_uxz, *_uyz, *_uzz
*/
) {
	int i, j, k, ii;
	int flag1 = 0, flag2 = 0;

	double _d, _p, _q;
	double _et, xi;
	double _u[12] = { 0.0 };
	double du[12] = { 0.0 };
	double _dua[12] = { 0.0 };
	double _dub[12] = { 0.0 };
	double _duc[12] = { 0.0 };

/**/
	if ( _z > 0.0 ) fprintf(stderr, "dc3d: Positive Z was given!\n");
/**/
	dccon0( _alpha, _dip );

/**/
	_d = _depth + _z;
	_p = _y * CD + _d * SD;
	_q = _y * SD - _d * CD;

	if ( (_x + _al1)*(_x - _al2) <= 0.0 ) flag1 = 1;
	if ( (_p + _aw1)*(_p - _aw2) <= 0.0 ) flag2 = 1;

/**/
	for ( k = 0; k < 2; k++ ) {
		if ( k == 0 ) _et = _p + _aw1;
		else _et = _p - _aw2;
	/**/
		for ( j = 0; j < 2; j++ ) {
			if ( j == 0 ) _xi = _x + _al1;
			else _xi = _x - _al2;
		/**/
			dccon2( &_xi, &_et, &_q, SD, CD );
		/**/
			if ( (flag1 && fabs(_q) < DBL_EPSILON && fabs(_et) < DBL_EPSILON) ||
				(flag2 && fabs(_q) < DBL_EPSILON && fabs(_xi) < DBL_EPSILON) )
				goto singular;
		/**/
			ua( _xi, _et, _q, _disl1, _disl2, _disl3, _dua );
		/**/
			ii = 0;
			for ( i = 0; i < 4; i++ ) {
				_du[ii] = -_dua[ii];
				ii++;
				_du[ii] = -_dua[ii]*CD + _dua[ii+1]*SD;
				ii++;
				_du[ii] = -_dua[ii-1]*SD - _dua[ii]*CD;
				ii++;
			}
		/**/
			_du[ 9] = -_du[9];
			_du[10] = -_du[10];
			_du[11] = -_du[11];
		/**/
			for ( i = 0; i < 12; i++ ) {
				if ( j + k != 1 )
					_u[i] += _du[i];
				if ( j + k == 1 )
					_u[i] -= _du[i];
			}
		}
	}

/**/
	_d = _depth - _z;
	_p = _y * CD + _d * SD;
	_q = _y * SD - _d * CD;
/**/
	for ( k = 0; k < 2; k++ ) {
		if ( k == 0 ) _et = _p + _aw1;
		else _et = _p - _aw2;
	/**/
		for ( j = 0; j < 2; j++ ) {
			if ( j == 0 ) _xi = _x + _al1;
			else _xi = _x - _al2;
		/**/
			dccon2( &_xi, &_et, &_q, SD, CD );
		/**/
			ua( _xi, _et, _q, _disl1, _disl2, _disl3, _dua );
			ub( _xi, _et, _q, _disl1, _disl2, _disl3, _dub );
			uc( _xi, _et, _q, _z, _disl1, _disl2, _disl3, _duc );
		/**/
			ii = 0;
			for ( i = 0; i < 4; i++ ) {
				_du[ii] = _dua[ii] + _dub[ii] + _z*_duc[ii];
				ii++;
				_du[ii] = (_dua[ii] + _dub[ii] + _z*_duc[ii])*CD - (_dua[ii+1] + _dub[ii+1] + _z*_duc[ii+1])*SD;
				ii++;
				_du[ii] = (_dua[ii-1] + _dub[ii-1] - _z*_duc[ii-1])*SD + (_dua[ii] + _dub[ii] - _z*_duc[ii])*CD;
				ii++;
			}
		/**/
			_du[ 9] += _duc[0];
			_du[10] += _duc[1]*CD - _duc[2]*SD;
			_du[11] -= _duc[1]*SD + _duc[2]*CD;
		/**/
			for ( i = 0; i < 12; i++ ) {
				if ( j + k != 1 )
					_u[i] += _du[i];
				if ( j + k == 1 )
					_u[i] -= _du[i];
			}
		}
	}

	memcpy(u, _u, sizeof(double)*12);
	return 0;

singular:
	memset(u, 0, sizeof(double)*12);
	return 1;
}


/**/
static void ua( double xi, double et, double q, double disl1, double disl2, double disl3, double u[12] )
{
	int i;

	double tmp1, tmp2, tmp3;
	double du[12];

	const double xy = xi * Y11;
	const double qx =  q * X11;
	const double qy =  q * Y11;
	const double half_tt = TT * 0.5;

	memset(u, 0, sizeof(double)*12);

/**/
	if ( fabs(disl1) > DBL_EPSILON ) {
	/**/
		tmp1  = ALP[1] * q;
		du[0] = half_tt      + tmp1 * xi * Y11;
		du[1] =                tmp1 / R;
		du[2] = ALP[0] * ALE - tmp1 * qy;
	/**/
		tmp2  = tmp1 * xi;
		du[3] = -ALP[0] * qy - tmp2 * xi * Y32;
		du[4] =              - tmp2 / R3;
		du[5] = ALP[0] * xy  + tmp2 * q * Y32;
	/**/
		tmp2   = X11 * 0.5;
		tmp3   = ALP[1] * xi;
		du[6]  = ALP[0] * xy * SD            + tmp3 * FY + D * tmp2;
		du[7]  =                               ALP[1] * EY;
		du[8]  = ALP[0] * (CD / R + qy * SD) - tmp1 * FY;
	/**/
		du[9]  = ALP[0] * xy * CD             + tmp3 * FZ + Y * tmp2;
		du[10] =                                ALP[1] * EZ;
		du[11] = -ALP[0] * (SD / R - qy * CD) - tmp1 * FZ;
	/**/
		tmp1 = disl1 / PI_D;
		for( i = 0; i < 12; i++ ) u[i] += tmp1 * du[i];
	}

/**/
	if ( fabs(disl2) > DBL_EPSILON ) {
	/**/
		tmp1  = ALP[1] * q;
		du[0] =                tmp1 / R;
		du[1] = half_tt      + tmp1 * et * X11;
		du[2] = ALP[0] * ALX - tmp1 * qx;
	/**/
		tmp2  = tmp1 / R3;
		du[3] =            - tmp2 * xi;
		du[4] = -qy * 0.5  - tmp2 * et;
		du[5] = ALP[0] / R + tmp2 * q;
	/**/
		tmp2  = ALP[0] * D * X11;
		tmp3  = ALP[0] * Y * X11;
		du[6] =                        ALP[1] * EY;
		du[7] = tmp2 + xy * 0.5 * SD + ALP[1] * et * GY;
		du[8] = tmp3                 - tmp1 * GY;
	/**/
		du[9]  =                        ALP[1] * EZ;
		du[10] = tmp3 + xy * 0.5 * CD + ALP[1] * et * GZ;
		du[11] = -tmp2                - tmp1 * GZ;
	/**/
		tmp1 = disl2 / PI_D;
		for( i = 0; i < 12; i++ ) u[i] += tmp1 * du[i];
	}

/**/
	if ( fabs(disl3) > DBL_EPSILON ) {
	/**/
		tmp1  = ALP[1] * q;
		du[0] = -ALP[0] * ALE - tmp1 * qy;
		du[1] = -ALP[0] * ALX - tmp1 * qx;
		du[2] = half_tt       - tmp1 * (et * X11 + xi * Y11);
	/**/
		tmp2  = tmp1 * q;
		du[3] = -ALP[0] * xy + tmp2 * xi * Y32;
		du[4] = -ALP[0] / R  + tmp2 / R3;
		du[5] = -ALP[0] * qy - tmp2 * q * Y32;
	/**/
		tmp2  = ALP[0] * Y * X11;
		tmp3  = ALP[0] * D * X11;
		du[6] = -ALP[0] * (CD / R + qy * SD) - tmp1 * FY;
		du[7] = -tmp2                        - tmp1 * GY;
		du[8] = tmp3 + ALP[0] * xy * SD      + tmp1 * HY;
	/**/
		du[9]  = ALP[0] * (SD / R - qy * CD) - tmp1 * FZ;
		du[10] = tmp3                        - tmp1 * GZ;
		du[11] = tmp2 + ALP[0] * xy * CD     + tmp1 * HZ;
	/**/
		tmp1 = disl3 / PI_D;
		for( i = 0; i < 12; i++ ) u[i] += tmp1 * du[i];
	}

	return;
}

/**/
static void ub ( double xi, double et, double q, double disl1, double disl2, double disl3, double u[12] )
{
	int i;

	double tmp;
	double du[12];

	const double xy = xi * Y11;
	const double qx =  q * X11;
	const double qy =  q * Y11;

	const double ir  = 1.0 / R;
	const double ir3 = ir * ir * ir;
	const double rd  = R + D;
	const double ird = 1.0 / rd;
	const double d11 = ir * ird;

	tmp = Y * ird * d11;
	const double aj2 = xi * tmp;
	const double aj5 = -(D * d11 + Y * tmp);

	double aj1, aj3, aj4, aj6;
	double ai1, ai2, ai3, ai4;
	double ak1, ak2, ak3, ak4;

/**/
	memset(u, 0, sizeof(double)*12);

/**/
	if ( fabs(CD) > DBL_EPSILON ) {
		if ( fabs(xi) < DBL_EPSILON ) {
			ai4 = 0.0;
		}
		else {
			tmp = sqrt(XI2 + Q2);
			ai4 = (xi * ird * SDCD + 2.0 * atan((et * (tmp + q * CD) + tmp * (R + tmp) * SD) / (xi * (R + tmp) * CD))) / CDCD;
		}
		tmp = 1.0 / CD;
		ai3 = (Y * CD * ird - ALE + SD * log(rd)) / CDCD;
		ak1 = xi * (d11 - Y11 * SD) * tmp;
		ak3 = (q * Y11 - Y * d11) * tmp;
		aj3 = (ak1 - aj2 * SD) * tmp;
		aj6 = (ak3 - aj5 * SD) * tmp;
	}
	else {
		tmp = ird * ird;
		ai3 = (et * ird + Y * q * tmp - ALE) * 0.5;
		ai4 = xi * Y * tmp * 0.5;
		ak1 = xi * q * ird * d11;
		ak3 = SD * ird * (XI2 * d11 - 1.0);
		aj3 = -xi * tmp * (Q2 * d11 - 0.5);
		aj6 = -Y * tmp * (XI2 * d11 - 0.5);
	}

	ai1 = -xi * ird * CD - ai4 * SD;
	ai2 = log(rd) + ai3 * SD;
	ak2 = ir + ak3 * SD;
	ak4 = xy * CD - ak1 * SD;
	aj1 = aj5 * CD - aj6 * SD;
	aj4 = -xy - aj2 * CD + aj3 * SD;

/**/
	if ( fabs(disl1) > DBL_EPSILON ) {
		tmp = ALP[2] * SD;
	/**/
		du[0] = -xi * qy - TT - tmp * ai1;
		du[1] = -q * ir       + tmp * Y * ird;
		du[2] = q * qy        - tmp * ai2;
	/**/
		du[3] = XI2 * q * Y32  - tmp * aj1;
		du[4] = xi * q * ir3   - tmp * aj2;
		du[5] = -xi * Q2 * Y32 - tmp * aj3;
	/**/
		du[6] = -xi * FY - D * X11 + tmp * (xy + aj4);
		du[7] = -EY                + tmp * (ir + aj5);
		du[8] = q * FY             - tmp * (qy + aj6);
	/**/
		du[9]  = -xi * FZ - Y * X11 + tmp * ak1;
		du[10] = -EZ                + tmp * Y * d11;
		du[11] = q * FZ             + tmp * ak2;
	/**/
		tmp = disl1 / PI_D;
		for ( i = 0; i < 12; i++ ) u[i] += tmp * du[i];
	}

/**/
	if ( fabs(disl2) > DBL_EPSILON ) {
		tmp = ALP[2] * SDCD;
	/**/
		du[0] = -q * ir       + tmp * ai3;
		du[1] = -et * qx - TT - tmp * xi * ird;
		du[2] = q * qx        + tmp * ai4;
	/**/
		du[3] = xi * q * ir3      + tmp * aj4;
		du[4] = et * q * ir3 + qy + tmp * aj5;
		du[5] = -Q2 * ir3         + tmp * aj6;
	/**/
		du[6] = -EY                + tmp * aj1;
		du[7] = -et * GY - xy * SD + tmp * aj2;
		du[8] = q * GY             + tmp * aj3;
	/**/
		du[9]  = -EZ                - tmp * ak3;
		du[10] = -et * GZ - xy * CD - tmp * xi * d11;
		du[11] = q * GZ             - tmp * ak4;
	/**/
		tmp = disl2 / PI_D;
		for ( i = 0; i < 12; i++ ) u[i] += tmp * du[i];
	}

/**/
	if ( fabs(disl3) > DBL_EPSILON ) {
		tmp = ALP[2] * SDSD;
	/**/
		du[0] = q * qy                 - tmp * ai3;
		du[1] = q * qx                 + tmp * xi * ird;
		du[2] = et * qx + xi * qy - TT - tmp * ai4;
	/**/
		du[3] = -xi * Q2 * Y32 - tmp * aj4;
		du[4] = -Q2 * ir3      - tmp * aj5;
		du[5] = q * Q2 * Y32   - tmp * aj6;
	/**/
		du[6] = q * FY  - tmp * aj1;
		du[7] = q * GY  - tmp * aj2;
		du[8] = -q * HY - tmp * aj3;
	/**/
		du[9]  = q * FZ  + tmp * ak3;
		du[10] = q * GZ  + tmp * xi * d11;
		du[11] = -q * HZ + tmp * ak4;
	/**/
		tmp = disl3 / PI_D;
		for ( i = 0; i < 12; i++ ) u[i] += tmp * du[i];
	}

	return;
}

/**/
static void uc ( double xi, double et, double q, double z, double disl1, double disl2, double disl3, double u[12] )
{
	int i;

	double tmp1, tmp2, tmp3;
	double du[12];

	const double xy = xi * Y11;
	/* const double qx =  q * X11; */
	const double qy =  q * Y11;

	const double ir   = 1.0 / R;
	const double ir2  = ir * ir;
	const double ir3  = ir * ir2;
	const double ir5  = ir2 * ir3;
	const double sdr3 = SD * ir3;
	const double cdr3 = CD * ir3;
	const double c    = D + z;
	const double qx53 = (8.0 * R2 + 9.0 * R * xi + 3.0 * XI2) * X11 * X11 * X11 * ir2 * q;

	tmp1 = q * CD - z;
	tmp2 = 3.0 * ir5;
	const double z32 = sdr3 - tmp1 * Y32;

	tmp3 = tmp2 * SD - tmp1 * (8.0 * R2 + 9.0 * R * et + 3.0 * ET2) * Y11 * Y11 * Y11 * ir2;
	const double y0  = Y11 - XI2 * Y32;
	const double z0  = z32 - XI2 * tmp3;

	tmp1 = q * Y32 * SD;
	const double xppy = (cdr3 + tmp1) * xi;
	const double xppz = (sdr3 - tmp1) * xi;
	const double qq  = z * Y32  + z32 + z0;

	tmp2 *= c;
	const double xqqy = (tmp2 * D - qq * SD) * xi;
	const double xqqz = (tmp2 * Y - qq * CD + q * Y32) * xi;
	const double cqr  = tmp2 * q;
	/* const double cqx = c * qx53; */
	tmp1 = c + D;
	/* cdr = (c + D) / R3 */
	const double cdrs = tmp1 * sdr3;
	const double cdrc = tmp1 * cdr3;
	const double yy0 = Y * ir3 - y0 * CD;

/**/
	memset(u, 0, sizeof(double)*12);

/**/
	if ( fabs(disl1) > DBL_EPSILON ) {
	/**/
		tmp1 = ALP[3] * CD;
		du[0] = tmp1 * xy                           - ALP[4] * xi * q * z32;
		du[1] = tmp1 * ir + ALP[3] * 2.0 * qy * SD) - ALP[4] * c * q * ir3;
		du[2] = tmp1 * qy                           - ALP[4] * (c * et * ir3 - z * Y11 + XI2 * z32);
	/**/
		du[3] = tmp1 * y0                                        - ALP[4] * q * z0;
		du[4] = -xi * (tmp1 * ir3 + ALP[3] * 2.0 * q * Y32 * SD) + ALP[4] * cqr * xi;
		du[5] = -tmp1 * xi * q * Y32                             + ALP[4] * xi * (3.0 * c * et * ir5 - qq);
	/**/
		du[6] = -tmp1 * xppy                                     - ALP[4] * xqqy;
		du[7] = ALP[3] * 2.0 * (D * sdr3 - y0 * SDSD) - Y * cdr3 - ALP[4] * (cdrs - et * ir3 - cqr * Y);
		du[8] = -ALP[3] * q * ir3 + yy0 * SD                     + ALP[4] * (cdrc + cqr * D - (y0 * CD + q * z0) * SD);
	/**/
		du[9]  = tmp1 * xppz                                      - ALP[4] * xqqz;
		du[10] = ALP[3] * 2.0 * (Y * sdr3 - y0 * SDCD) + D * cdr3 - ALP[4] * (cdrc + cqr * D);
		du[11] = yy0 * CD                                         - ALP[4] * (cdrs - cqr * Y - y0 * SDSD + q * z0 * CD);
	/**/
		tmp = disl1 / PI_D;
		for( i = 0; i < 12; i++ ) u[i] += tmp * du[i];
	}

/**/
	if ( fabs(disl2) > DBL_EPSILON ) {
	/**/
		du[0] = ALP[3] * CD * ir - qy * SD - ALP[4] * c * q * ir3;
		du[1] = ALP[3] * Y * X11           - ALP[4] * c * et * q * X32;
		du[2] = -D * X11 - xy * SD         - ALP[4] * c * (X11 - Q2 * X32);
	/**/
		du[3] = -ALP[3] * xi * cdr3 + ALP[4] * cqr * xi + xi * q * Y32 * SD;
		du[4] = -ALP[3] * Y * ir3   + ALP[4] * cqr * et;
		du[5] = D * ir3 - y0 * SD   + ALP[4] * c * ir3 * (1.0 - 3.0 * Q2 * ir2);
	/**/
		du[6] = -ALP[3] * et * ir3 + y0 * SDSD - ALP[4] * (cdrs - cqr * Y);
		du[7] = ALP[3] * (X11 - Y * Y * X32)   - ALP[4] * c * ((D + 2.0 * q * CD) * X32 - Y * et * qx53);
		du[8] = xppy * SD + Y * D * X32        + ALP[4] * c * ((Y + 2.0 * q * SD) * X32 - Y * q * qx53);
	/**/
		du[9]  = -q * ir3 + y0 * SDCD           - ALP[4] * (cdrc + cqr * D);
		du[10] = ALP[3] * Y * D * X32           - ALP[4] * c * ((Y - 2.0 * q * SD) * X32 + D * et * qx53);
		du[11] = -xppz * SD + X11 - D * D * X32 - ALP[4] * c * ((D - 2.0 * q * CD) * X32 - D * q * qx53);
	/**/
		tmp = disl2 / PI_D;
		for( i = 0; i < 12; i++ ) u[i] += tmp * du[i];
	}

/**/
	if ( fabs(disl3) > DBL_EPSILON ) {
		tmp = ALP[4] * c;
	/**/
		du[0] = -ALP[3] * (SD * ir + qy * CD)    - ALP[4] * (z * Y11 - Q2 * z32);
		du[1] = ALP[3] * 2.0 * xy * SD + D * X11 - ALP[4] * c * (X11 - Q2 * X32);
		du[2] = ALP[3] * (Y * X11 + xy * CD)     + ALP[4] * q * (c * et * X32 + xi * z32);
	/**/
		du[3] = ALP[3] * xi * ir3 * SD + xi * q * Y32 * CD + ALP[4] * xi * (3.0 * c * et * ir5 - 2.0 * z32 - z0);
		du[4] = ALP[3] * 2.0 * y0 * SD - D * ir3           + ALP[4] * c * ir3 * (1.0 - 3.0 * Q2 * ir2);
		du[5] = -ALP[3] * yy0                              - ALP[4] * (cqr * et - q * z0);
	/**/
		du[6] = ALP[3] * (q * ir3 + y0 * SDCD)            + ALP[4] * (z * ir3 * CD + cqr * D - q * z0 * SD);
		du[7] = -ALP[3] * 2.0 * xppy * SD - Y * D * X32   + ALP[4] * c * ((Y + 2.0 * q * SD) * X32 - Y * q * qx53);
		du[8] = -ALP[3] * (xppy * CD - X11 + Y * Y * X32) + ALP[4] * (c * ((D + 2.0 * q * CD) * X32 - Y * et * qx53) + xqqy);
	/**/
		du[9]  = -et * ir3 + y0 * CDCD                        - ALP[4] * (z * ir3 * SD - cqr * Y - y0 * SDSD + q * z0 * CD);
		du[10] = ALP[3] * 2.0 * xppz * SD - X11 + D * D * X32 - ALP[4] * c * ((D - 2.0 * q * CD) * X32 - D * q * qx53);
		du[11] = ALP[3] * (xppz * CD + Y * D * X32)           + ALP[4] * (c * ((Y - 2.0 * q * SD) * X32 + D * et * qx53) + xqqz);
	/**/
		tmp = disl3 / PI_D;
		for( i = 0; i < 12; i++ ) u[i] += tmp * du[i];
	}

	return;
}

/**/
static void dccon0 ( double alpha, double dip )
{
	ALP[1] = alpha * 0.5;
	ALP[0] = 0.5 - ALP[1];   /* (1.0 - _alpha) * 0.5 */
	ALP[3] = 1.0 - alpha;
	ALP[2] = ALP[3] / alpha; /* (1.0 - _alpha) / _alpha */
	ALP[4] = alpha;

	dip *= DEG2RAD;
	SD = sin(dip);
	CD = cos(dip);

	if ( fabs(CD) < L_EPSILON ) {
		CD = 0.0;
		SD = SD > 0.0 ? 1.0 : -1.0;
	}

	SDSD = SD * SD;
	CDCD = CD * CD;
	SDCD = SD * CD;
	S2D  = SDCD * 2.0;
	C2D  = CDCD - SDSD;

	return;
}

/**/
static void dccon2 ( double *xi, double *et, double *q, double sd, double cd )
{
	double tmp1, tmp2, tmp3;

	*xi = fabs(*xi) < L_EPSILON ? 0.0 : *xi;
	*et = fabs(*et) < L_EPSILON ? 0.0 : *et;
	*q  = fabs(*q)  < L_EPSILON ? 0.0 : *q;

	XI2 = *xi * *xi;
	ET2 = *et * *et;
	Q2  = *q * *q;
	R2  = XI2 + ET2 + Q2;
	R   = sqrt(R2);
	if ( R < DBL_EPSILON ) return;
	R3  = R * R2;
	R5  = R3 * R2;
	Y   = *et * cd + *q * sd;
	D   = *et * sd - *q * cd;
	TT  = fabs(*q) < DBL_EPSILON ? 0.0 : atan(*xi * *et / (*q * R));

	if ( *xi < 0.0 && fabs(*q) < DBL_EPSILON && fabs(*et) < DBL_EPSILON ) {
		ALX = -log(R - *xi);
		X11 = 0.0;
		X32 = 0.0;
	}
	else {
		tmp1 = R + *xi;
		ALX = log(tmp1);
		X11 = 1.0 / (R * tmp1);
		X32 = (R + tmp1) * X11 * X11 / R;
	}

	if ( *et < 0.0 && fabs(*q) < DBL_EPSILON && fabs(*xi) < DBL_EPSILON ) {
		ALE = -log(R - *et);
		Y11 = 0.0;
		Y32 = 0.0;
	}
	else {
		tmp1 = R + *et;
		ALE = log(tmp1);
		Y11 = 1.0 / (R * tmp1);
		Y32 = (R + tmp1) * Y11 * Y11 / R;
	}

	tmp1 = Y / R3;
	tmp2 = D / R3;
	EY = sd / R - *q * tmp1;
	EZ = cd / R + *q * tmp2;
	tmp3 = XI2 * Y32;
	FY = tmp2 + tmp3 * sd;
	FZ = tmp1 + tmp3 * cd;

	tmp1 = X11 * 2.0;
	tmp2 = Y * *q * X32;
	tmp3 = D * *q * X32;
	GY = tmp1 * sd - tmp2;
	GZ = tmp1 * cd + tmp3;
	tmp1 = *xi * *q * Y32;
	HZ = tmp2 + tmp1 * cd;
	HY = tmp3 + tmp1 * sd;

	return;
}
