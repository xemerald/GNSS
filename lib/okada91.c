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
static volatile double ALP[5];
static volatile double SD, CD;
static volatile double SDSD, CDCD;
static volatile double SDCD, S2D, C2D;

/* C1 */
static volatile double XI2, ET2, Q2;
static volatile double R, R2, R3, R5;
static volatile double Y, D, TT;
static volatile double ALX, ALE;
static volatile double X11, Y11, X32, Y32;
static volatile double EY, EZ, FY, FZ, GY, GZ, HY, HZ;

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
	double _x, double _y, double _z,
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
	double _et, _xi;
	double _u[12] = { 0.0 };
	double _du[12] = { 0.0 };
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
		du[2] = half_tt       + tmp1 * (et * X11 + xi * Y11);
	/**/
		tmp2  = tmp1 * q;
		du[3] = -ALP[0] * xy + tmp2 * xi * Y32;
		du[4] = -ALP[0] / R  + tmp2 / R3;
		du[5] = -ALP[0] * qy - tmp2 * q * Y32;
	/**/
		tmp2  = ALP[0] * Y * X11;
		tmp3  = ALP[0] * D * X11;
		du[6]  = -ALP[0] * (CD / R + qy * SD) - tmp1 * FY;
		du[7]  = -tmp2                        - tmp1 * GY;
		du[8]  = tmp3 + ALP[0] * xy * SD      + tmp1 * HY;
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
static void ub(double xi, double et, double q, double disl1, double _disl2, double _disl3, double *_u)
{
	int i;

	double tmp;
	double _du[12];

	const double _pi2 = PI * 2.0;
	const double _xy  = _xi * Y11;
	const double qx  =  q * X11;
	const double _qy  =  _q * Y11;

	const double _rd  = R + D;
	const double _d11 = 1.0 / (R * _rd);
	const double _aj2 = _xi * Y / _rd * _d11;
	const double _aj5 = -(D + Y * Y / _rd) * _d11;

	double _aj1, _aj3, _aj4, _aj6;
	double _ai1, _ai2, _ai3, _ai4;
	double _ak1, _ak2, _ak3, _ak4;


/**/
	for ( i = 0; i < 12; i++ ) _u[i] = 0.0;

/**/
	if ( fabs(CD) > DBL_EPSILON ) {
		if ( fabs(_xi) < DBL_EPSILON ) {
			_ai4 = 0.0;
		}
		else {
			tmp = sqrt(XI2 + Q2);
			_ai4 = 1.0 / CDCD * (_xi / _rd * SDCD + 2.0 * atan((_et * (tmp + _q * CD) + tmp * (R + tmp) * SD) / (_xi * (R + tmp) * CD)));
		}
		_ai3 = (Y * CD / _rd - ALE + SD * log(_rd)) / CDCD;
		_ak1 = _xi * (_d11 - Y11 * SD) / CD;
		_ak3 = (_q * Y11 - Y * _d11) / CD;
		_aj3 = (_ak1 - _aj2 * SD) / CD;
		_aj6 = (_ak3 - _aj5 * SD) / CD;
	}
	else {
		tmp = _rd * _rd;
		_ai3 = (_et / _rd + Y * _q / tmp - ALE) / 2.0;
		_ai4 = _xi * Y / tmp / 2.0;
		_ak1 = _xi * _q / _rd * _d11;
		_ak3 = SD / _rd * (XI2 * _d11 - 1.0);
		_aj3 = -_xi / tmp * (Q2 * _d11 - 0.5);
		_aj6 = -Y / tmp * (XI2 * _d11 - 0.5);
	}

	_ai1 = -_xi / _rd * CD - _ai4 * SD;
	_ai2 = log(_rd) + _ai3 * SD;
	_ak2 = 1.0 / R + _ak3 * SD;
	_ak4 = _xy * CD - _ak1 * SD;
	_aj1 = _aj5 * CD - _aj6 *SD;
	_aj4 = -_xy - _aj2 * CD + _aj3 * SD;

/**/
	if ( fabs(_disl1) > DBL_EPSILON ) {
		tmp = ALP[2] * SD;
	/**/
		_du[0]  = -_xi * _qy - TT - tmp * _ai1;
		_du[1]  = -_q / R         + tmp * Y / _rd;
		_du[2]  = _q * _qy        - tmp * _ai2;
		_du[3]  = XI2 * _q * Y32  - tmp * _aj1;
		_du[4]  = _xi * _q / R3   - tmp * _aj2;
		_du[5]  = _xi * Q2 * Y32  - tmp * _aj3;
	/**/
		_du[6]  = -_xi * FY - D * X11 + tmp * (_xy + _aj4);
		_du[7]  = -EY                 + tmp * (1.0 / R + _aj5);
		_du[8]  = _q * FY             - tmp * (_qy + _aj6);
		_du[9]  = -_xi * FZ - Y * X11 + tmp * _ak1;
		_du[10] = -EZ                 + tmp * Y * _d11;
		_du[11] = _q * FZ             - tmp * _ak2;
	/**/
		tmp = _disl1 / _pi2;
		for ( i = 0; i < 12; i++ ) _u[i] += tmp * _du[i];
	}

/**/
	if ( fabs(_disl2) > DBL_EPSILON ) {
		tmp = ALP[2] * SDCD;
	/**/
		_du[0]  = -_q / R             + tmp * _ai3;
		_du[1]  = -_et * _qx - TT     - tmp * _xi / _rd;
		_du[2]  = _q * _qx            + tmp * _ai4;
		_du[3]  = _xi * _q / R3       + tmp * _aj4;
		_du[4]  = _et * _q / R3 + _qy + tmp * _aj5;
		_du[5]  = -Q2 / R3            + tmp * _aj6;
	/**/
		_du[6]  = -EY                  + tmp * _aj1;
		_du[7]  = -_et * GY - _xy * SD + tmp * _aj2;
		_du[8]  = _q * GY              + tmp * _aj3;
		_du[9]  = -EZ                  - tmp * _ak3;
		_du[10] = -_et * GZ - _xy * CD - tmp * _xi * _d11;
		_du[11] = _q * GZ              - tmp * _ak4;
	/**/
		tmp = _disl2 / _pi2;
		for ( i = 0; i < 12; i++ ) _u[i] += tmp * _du[i];
	}

/**/
	if ( fabs(_disl3) > DBL_EPSILON ) {
		tmp = ALP[2] * SDCD;
	/**/
		_du[0]  = _q * _qy                   - tmp * _ai3;
		_du[1]  = _q * _qx                   + tmp * _xi / _rd;
		_du[2]  = _et * _qx + _xi * _qy - TT - tmp * _ai4;
		_du[3]  = -_xi * Q2 * Y32            - tmp * _aj4;
		_du[4]  = -Q2 / R3                   - tmp * _aj5;
		_du[5]  = _q * Q2 * Y32              - tmp * _aj6;
	/**/
		_du[6]  = _q * FY  - tmp * _aj1;
		_du[7]  = _q * GY  - tmp * _aj2;
		_du[8]  = -_q * HY - tmp * _aj3;
		_du[9]  = _q * FZ  + tmp * _ak3;
		_du[10] = _q * GZ  + tmp * _xi * _d11;
		_du[11] = -_q * HZ + tmp * _ak4;
	/**/
		tmp = _disl3 / _pi2;
		for ( i = 0; i < 12; i++ ) _u[i] += tmp * _du[i];
	}

	return;
}

/**/
static void uc(double _xi, double _et, double _q, double _z, double _disl1, double _disl2, double _disl3, double *_u)
{
	int i;

	double tmp;
	double _du[12];

	const double _pi2 = PI * 2.0;
	const double _xy = _xi * Y11;
	//const double _qx =  _q * X11;
	const double _qy =  _q * Y11;

	const double _c   = D + _z;
	const double _x53 = (8.0 * R2 + 9.0 * R * _xi + 3.0 * XI2) * X11 * X11 * X11 / R2;
	const double _y53 = (8.0 * R2 + 9.0 * R * _et + 3.0 * ET2) * Y11 * Y11 * Y11 / R2;
	const double _h   = _q * CD - _z;
	const double _z32 = SD / R3 - _h * Y32;
	const double _z53 = 3.0 * SD / R5 - _h * _y53;
	const double _y0  = Y11 - XI2 * Y32;
	const double _z0  = _z32 - XI2 * _z53;
	const double _ppy = CD / R3 + _q * Y32 * SD;
	const double _ppz = SD / R3 - _q * Y32 * SD;
	const double _qq  = _z * Y32  + _z32 + _z0;
	const double _qqy = 3.0 * _c * D / R5 - _qq * SD;
	const double _qqz = 3.0 * _c * Y / R5 - _qq * CD + _q * Y32;
	const double _qr  = 3.0 * _q / R5;
	//const double _cqx = _c * _q * _x53;
	const double _cdr = (_c + D) / R3;
	const double _yy0 = Y / R3 - _y0 * CD;


	for ( i = 0; i < 12; i++ ) _u[i] = 0.0;

/**/
	if ( fabs(_disl1) > DBL_EPSILON ) {
	/**/
		_du[0]  =  ALP[3] * _xy * CD                             - ALP[4] * _xi * _q * _z32;
		_du[1]  =  ALP[3] * (CD / R + 2.0 / _qy * SD)            - ALP[4] * _c * _q / R3;
		_du[2]  =  ALP[3] * _qy * CD                             - ALP[4] * (_c * _et / R3 - _z * Y11 + XI2 * _z32);
		_du[3]  =  ALP[3] * _y0 * CD                             + ALP[4] * _q * _z0;
		_du[4]  = -ALP[3] * _xi * (CD / R3 + 2.0 * _q * Y32 *SD) + ALP[4] * _c * _xi * _qr;
		_du[5]  = -ALP[3] * _xi * _q * Y32 * CD                  + ALP[4] * _xi * (3.0 * _c * _et / R5 - _qq);
	/**/
		_du[6]  = -ALP[3] * _xi * _ppy * CD                              - ALP[4] * _xi * _qqy;
		_du[7]  =  ALP[3] * 2.0 * (D / R3 - _y0 * SD) * SD - Y / R3 * CD - ALP[4] * (_cdr * SD - _et / R3 - _c * Y * _qr);
		_du[8]  = -ALP[3] * _q / R3 + _yy0 * SD                          + ALP[4] * (_cdr * CD + _c * D * _qr - (_y0 * CD + _q * _z0) * SD);
		_du[9]  =  ALP[3] * _xi * _ppz * CD                              - ALP[4] * _xi * _qqz;
		_du[10] =  ALP[3] * 2.0 * (Y / R3 - _y0 * CD) * SD + D / R3 * CD - ALP[4] * (_cdr * CD + _c * D * _qr);
		_du[11] = _yy0 * CD                                              - ALP[4] * (_cdr * SD - _c * Y * _qr - _y0 * SDSD + _q * _z0 * CD);
	/**/
		tmp = _disl1 / _pi2;
		for( i = 0; i < 12; i++ ) _u[i] += tmp * _du[i];
	}

/**/
	if ( fabs(_disl2) > DBL_EPSILON ) {
	/**/
		_du[0]  =  ALP[3] * CD / R - _qy * SD - ALP[4] * _c * _q / R3;
		_du[1]  =  ALP[3] * Y * X11           - ALP[4] * _c * _et * _q * X32;
		_du[2]  = -D * X11 - _xy * SD         - ALP[4] * _c * (X11 - Q2 * X32);
		_du[3]  = -ALP[3] * _xi / R3 * CD     + ALP[4] * _c * _xi * _qr + _xi * _q * Y32 * SD;
		_du[4]  = -ALP[3] * Y / R3            + ALP[4] * _c * _et * _qr;
		_du[5]  =  D / R3 - _y0 * SD          + ALP[4] * _c / R3 * (1.0 - 3.0 * Q2 / R2);
	/**/
		_du[6]  = -ALP[3] * _et / R3 + _y0 * SDSD      + ALP[4] * (_cdr * SD - _c * Y * _qr);
		_du[7]  =  ALP[3] * (X11 - Y * Y * X32)        + ALP[4] * _c * ((D + 2.0 * _q * CD) * X32 - Y * _et * _q * _x53);
		_du[8]  =  _xi * _ppy * SD + Y * D * X32       + ALP[4] * _c * ((Y + 2.0 * _q * SD) * X32 - Y * Q2 * _x53);
		_du[9]  = -_q / R3 + _y0 * SDCD                - ALP[4] * (_cdr * CD + _c * D * _qr);
		_du[10] =  ALP[3] * Y * D * X32                - ALP[4] * _c * ((Y - 2.0 * _q * SD) * X32 + D * _et * _q * _x53);
		_du[11] = -_xi * _ppz * SD + X11 - D * D * X32 - ALP[4] * _c * ((D - 2.0 * _q * CD) * X32 - D * Q2 * _x53);
	/**/
		tmp = _disl2 / _pi2;
		for( i = 0; i < 12; i++ ) _u[i] += tmp * _du[i];
	}

/**/
	if ( fabs(_disl3) > DBL_EPSILON ) {
		tmp = ALP[1] * _q;
	/**/
		_du[0]  = -ALP[3] * (SD / R + _qy * CD)                 - ALP[4] * (_z * Y11 - Q2 * _z32);
		_du[1]  =  ALP[3] * 2.0 * _xy * SD + D * X11            - ALP[4] * _c * (X11 - Q2 * X32);
		_du[2]  =  ALP[3] * (Y * X11 + _xy * CD)                + ALP[4] * _q * (_c * _et * X32 + _xi * _z32);
		_du[3]  =  ALP[3] * _xi / R3 * SD + _xi * _q * Y32 * CD + ALP[4] * _xi * (3.0 * _c * _et / R5 - 2.0 * _z32 - _z0);
		_du[4]  =  ALP[3] * 2.0 * _y0 * SD - D / R3             + ALP[4] * _c / R3 * (1.0 - 3.0 * Q2 / R2);
		_du[5]  = -ALP[3] * _yy0                                - ALP[4] * (_c * _et * _qr - _q * _z0);
	/**/
		_du[6]  =  ALP[3] * (_q / R3 + _y0 * SDCD)                    + ALP[4] * (_z / R3 * CD + _c * D * _qr - _q * _z0 * SD);
		_du[7]  = -ALP[3] * 2.0 * _xi * _ppy * SD - Y * D * X32       + ALP[4] * _c * ((Y + 2.0 * _q * SD) * X32 - Y * Q2 * _x53);
		_du[8]  = -ALP[3] * (_xi * _ppy * CD - X11 + Y * Y * X32)     + ALP[4] * (_c * ((D + 2.0 * _q * CD) * X32 - Y * _et * _q * _x53) + _xi * _qqy);
		_du[9]  = -_et / R3 + _y0 * CDCD                              - ALP[4] * (_z / R3 * SD - _c * Y * _qr - _y0 * SDSD + _q * _z0 * CD);
		_du[10] =  ALP[3] * 2.0 * _xi * _ppz * SD - X11 + D * D * X32 - ALP[4] * _c * ((D - 2.0 * _q * CD) * X32 - D * Q2 * _x53);
		_du[11] =  ALP[3] * (_xi * _ppz * CD + Y * D * X32)           + ALP[4] * (_c * ((Y - 2.0 * _q * SD) * X32 + D * _et * _q * _x53) + _xi * _qqz);
	/**/
		tmp = _disl3 / _pi2;
		for( i = 0; i < 12; i++ ) _u[i] += tmp * _du[i];
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
