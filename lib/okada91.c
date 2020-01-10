#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "../include/constant.h"

/**/
static void dccon0(double, double);
static void dccon2(double, double, double, double, double);
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
int dc3d(double alpha, double x, double y, double z, double depth, double dip,
	double al1, double al2, double aw1, double aw2, double disl1, double disl2, double disl3,
	double ux, double uy, double uz, double uxx, double uyx, double uzx,
	double uxy,	double uyy, double uzy, double uxz, double uyz, double uzz)
{
	int i;
	double u[12], du[12];
	double dua[12], dub[12], duc[12];

	if ( z > 0.0 ) fprintf(stderr, "dc3d: Positive Z was given!\n");

	i = 12 * sizeof(double);
	memset(u,   0, i);
	memset(du,  0, i);
	memset(dua, 0, i);
	memset(dub, 0, i);
	memset(duc, 0, i);

	dccon0(alpha, dip);

	return 0;
}


/**/
static void ua(double _xi, double _et, double _q, double _disl1, double _disl2, double _disl3, double *_u)
{
	int i;

	double tmp;
	double _du[12];

	const double _pi2 = PI * 2.0;
	const double _xy = _xi * Y11;
	const double _qx =  _q * X11;
	const double _qy =  _q * Y11;

	for ( i = 0; i < 12; i++ ) _u[i] = 0.0;

/**/
	if ( fabs(_disl1) > DBL_EPSILON ) {
	/**/
		_du[0]  =  TT / 2.0     + ALP[1] * _xi * _qy;
		_du[1]  =                 ALP[1] * _q / R;
		_du[2]  =  ALP[0] * ALE - ALP[1] * _q * _qy;
		_du[3]  = -ALP[0] * _qy - ALP[1] * XI2 * _q * Y32;
		_du[4]  =               - ALP[1] * _xi * _q / R3;
		_du[5]  =  ALP[0] * _xy + ALP[1] * _xi * Q2 * Y32;
	/**/
		_du[6]  =  ALP[0] * _xy * SD            + ALP[1] * _xi * FY + D / 2.0 * X11;
		_du[7]  =                                 ALP[1] * EY;
		_du[8]  =  ALP[0] * (CD / R + _qy * SD) - ALP[1] * _q * FY;
		_du[9]  =  ALP[0] * _xy * CD            + ALP[1] * _xi * FZ + Y / 2.0 * X11;
		_du[10] =                                 ALP[1] * EZ;
		_du[11] = -ALP[0] * (SD / R - _qy * CD) - ALP[1] * _q * FZ;
	/**/
		tmp = _disl1 / _pi2;
		for( i = 0; i < 12; i++ ) _u[i] += tmp * _du[i];
	}

/**/
	if ( fabs(_disl2) > DBL_EPSILON ) {
	/**/
		_du[0]  =                 ALP[1] * _q / R;
		_du[1]  =  TT / 2.0     + ALP[1] * _et * _qx;
		_du[2]  =  ALP[0] * ALX - ALP[1] * _q * _qx;
		_du[3]  =               - ALP[1] * _xi * _q / R3;
		_du[4]  = -_qy / 2.0    - ALP[1] * _et * _q / R3;
		_du[5]  =  ALP[0] / R   + ALP[1] * Q2 / R3;
	/**/
		_du[6]  =                                      ALP[1] * EY;
		_du[7]  =  ALP[0] * D * X11 + _xy / 2.0 * SD + ALP[1] * _et * GY;
		_du[8]  =  ALP[0] * Y * X11                  - ALP[1] * _q * GY;
		_du[9]  =                                      ALP[1] * EZ;
		_du[10] =  ALP[0] * Y * X11 + _xy / 2.0 * CD + ALP[1] * _et * GZ;
		_du[11] = -ALP[0] * D * X11                  - ALP[1] * _q * GZ;
	/**/
		tmp = _disl2 / _pi2;
		for( i = 0; i < 12; i++ ) _u[i] += tmp * _du[i];
	}

/**/
	if ( fabs(_disl3) > DBL_EPSILON ) {
		tmp = ALP[1] * _q;
	/**/
		_du[0]  = -ALP[0] * ALE - tmp * _qy;
		_du[1]  = -ALP[0] * ALX - tmp * _qx;
		_du[2]  =  TT / 2.0     + ALP[1] * (_et * _qx + _xi * _qy);
		_du[3]  = -ALP[0] * _xy + ALP[1] * _xi * Q2 * Y32;
		_du[4]  = -ALP[0] / R   + ALP[1] * Q2 / R3;
		_du[5]  = -ALP[0] * _qy - tmp * Q2 * Y32;
	/**/
		_du[6]  = -ALP[0] * (CD / R + _qy * SD)  - tmp * FY;
		_du[7]  = -ALP[0] * Y * X11              - tmp * GY;
		_du[8]  =  ALP[0] * (D * X11 + _xy * SD) + tmp * HY;
		_du[9]  =  ALP[0] * (SD / R - _qy * CD)  - tmp * FZ;
		_du[10] =  ALP[0] * D * X11              - tmp * GZ;
		_du[11] =  ALP[0] * (Y * X11 + _xy * CD) + tmp * HZ;
	/**/
		tmp = _disl3 / _pi2;
		for( i = 0; i < 12; i++ ) _u[i] += tmp * _du[i];
	}

	return;
}

/**/
static void ub(double _xi, double _et, double _q, double _disl1, double _disl2, double _disl3, double *_u)
{
	int i;

	double tmp;
	double _du[12];

	const double _pi2 = PI * 2.0;
	const double _xy  = _xi * Y11;
	const double _qx  =  _q * X11;
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
		tmp = _disl2 / _pi2
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
		_du[9]  = _q*FZ    + tmp * _ak3;
		_du[10] = _q*GZ;   + tmp * _xi * _d11;
		_du[11] = -_q*HZ   + tmp * _ak4;
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
	const double _qx =  _q * X11;
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
	const double _cqx = _c * _q * _x53;
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
static void dccon0(double alpha, double dip)
{
	ALP[0] = (1.0 - alpha) * 0.5;
	ALP[1] = alpha * 0.5;
	ALP[2] = (1.0 - alpha) / alpha;
	ALP[3] = 1.0 - alpha;
	ALP[4] = alpha;

	dip = dip * (PI / 180.0);
	SD = sin(dip);
	CD = cos(dip);

	if ( fabs(CD) < DBL_EPSILON ) {
		CD = 0.0;
		if ( SD > 0.0 ) SD = 1.0;
		else SD = -1.0;
	}

	SDSD = SD * SD;
	CDCD = CD * CD;
	SDCD = SD * CD;
	S2D  = SDCD * 2.0;
	C2D  = CDCD - SDSD;

	return;
}

/**/
static void dccon2(double xi, double et, double q, double sd, double cd)
{
	double tmp;

	xi = xi < DBL_EPSILON ? 0.0 : xi;
	et = et < DBL_EPSILON ? 0.0 : et;
	q  = q  < DBL_EPSILON ? 0.0 : q;

	XI2 = xi * xi;
	ET2 = et * et;
	Q2  = q * q;
	R2  = XI2 + ET2 + Q2;
	R   = sqrt(R2);
	if ( R < DBL_EPSILON ) return;
	R3  = R * R2;
	R5  = R3 * R2;
	Y   = et * cd + q * sd;
	D   = er * sd - q * cd;

	if ( fabs(q) < DBL_EPSILON )
		TT = 0.0;
	else
		TT = atan(xi*et / (q*R));

	if ( xi < 0.0 && fabs(q) < DBL_EPSILON && fabs(et) < DBL_EPSILON ) {
		ALX = -log(R - xi);
		X11 = 0.0;
		X32 = 0.0;
	}
	else {
		tmp = R + xi;
		ALX = log(tmp);
		X11 = 1.0 / (R * tmp);
		X32 = (R + tmp)*X11*X11 / R;
	}

	if ( et < 0.0 && fabs(q) < DBL_EPSILON && fabs(xi) < DBL_EPSILON ) {
		ALE = -log(R - et);
		Y11 = 0.0;
		Y32 = 0.0;
	}
	else {
		tmp = R + et;
		ALE = log(tmp);
		Y11 = 1.0 / (R * tmp);
		Y32 = (R + tmp)*Y11*Y11 / R;
	}

	EY = sd / R - Y * q / R3;
	EZ = cd / R + D * q / R3;
	FY = D / R3 + XI2 * Y32 * sd;
	FZ = Y / R3 + XI2 * Y32 * cd;
	GY = 2.0 * X11 * sd - Y * q * X32;
	GZ = 2.0 * X11 * cd - D * q * X32;
	HY = D * q * X32 + xi * q * Y32 * sd;
	HZ = Y * q * X32 + xi * q * Y32 * cd;

	return;
}
