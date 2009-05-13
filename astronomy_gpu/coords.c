#include <math.h>

#include "coords.h"
#include "pi_constants.h"

double d_get_incl(int wedge) {
	double incl = (wedge * stripe_separation) - 57.5 + survey_center_dec;
	if (wedge > 46) incl -= 180.0;
	else incl -= 0.0;
	return incl; 
}

double d_get_incl_rad(int wedge) {
	return d_get_incl(wedge) * D_DEG2RAD;
}

float f_get_incl(int wedge) {
	float incl = (wedge * stripe_separation) - 57.5 + survey_center_dec;
	if (wedge > 46) incl -= 180.0;
	else incl -= 0.0;
	return incl;
}

float f_get_incl_rad(int wedge) {
	return f_get_incl(wedge) * F_DEG2RAD;
}

void d_lbr2xyz(const double* lbr, double* xyz) {
	double r0, sinb, sinl, cosb, cosl, zp, d;

	r0 = 8.5;
	sinb = sin(lbr[1]);
	sinl = sin(lbr[0]);
	cosb = cos(lbr[1]);
	cosl = cos(lbr[0]);

	xyz[2] = lbr[2] * sinb;
	zp = lbr[2] * cosb;
	d = sqrt( r0 * r0 + zp * zp - 2 * r0 * zp * cosl);
	xyz[0] = (zp * zp - r0 * r0 - d * d) / (2 * r0);
	xyz[1] = zp * sinl;
}
