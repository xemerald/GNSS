#pragma once

typedef struct {
	double length;
	double width;
	double depth;
	double dip;
	double strike;
	double east;
	double north;
} FAULT_GEOM;

FAULT_GEOM transfault_bottom( FAULT_GEOM );
int fault2patch( const FAULT_GEOM, const int, const int, FAULT_GEOM ** );
