#pragma once

typedef struct {
	double length;
	double width;
	double depth;
	double dip;
	double strike;
	double east;
	double north;
	double disl1;
	double disl2;
	double disl3;
} FAULT_MODEL;

typedef struct {
	double x;
	double y;
	double z;
} STA_COOR;

/* */
int disloc3d ( FAULT_MODEL *, STA_COOR *, double, double, double [], double [], double [] );
/* */
int dc3d ( double, double, double, double, double, double, double, double,
	double, double, double, double, double, double [] );
