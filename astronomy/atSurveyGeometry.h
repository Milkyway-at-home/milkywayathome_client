#ifndef ATSURVEYGEOMETRY_H
#define ATSURVEYGEOMETRY_H

#define D2PI 6.2831853071795864769252867665590057683943387987502
#define DPI 3.1415926535897932384626433832795028841971693993751
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
                                        :(A)+(B)*floor(-(A)/(B))):(A))
#define dsign(A,B) ((B)<0.0?-(A):(A))

#define at_stripeSeparation  2.5;
#define at_surveyCenterRa 185.0;
#define at_surveyCenterDec 32.5;
#define at_deg2Rad DPI/180.0;
#define at_rad2Deg 180.0/DPI;

void atGCToEq (
	       double amu,	/* IN */
	       double anu,	/* IN */
	       double *ra,	/* OUT */
	       double *dec,	/* OUT */
	       double anode,	/* IN */
	       double ainc	/* IN */
	       );

void atEqToGal (
		double ra,	/* IN */
		double dec,	/* IN */
		double *glong,	/* OUT: Galactic longitude */
		double *glat	/* OUT: Galactic latitude */
		);

void atBound (
	      double *angle,	/* MODIFIED -- the angle to bound */
	      double min,	/* IN -- inclusive minimum value */
	      double max	/* IN -- exclusive maximum value */
	      );

void atBound2(
	      double *theta,	/* MODIFIED -- the -90 to 90 angle */
	      double *phi	/* MODIFIED -- the 0 to 360 angle */
	      );

void slaDcc2s ( double v[3], double *a, double *b );

void slaDimxv ( double dm[3][3], double va[3], double vb[3] );

void slaDcs2c ( double a, double b, double v[3] );

void slaDmxv ( double dm[3][3], double va[3], double vb[3] );

double slaDrange ( double angle );

double slaDranrm ( double angle );

void slaEqgal ( double dr, double dd, double *dl, double *db );

#endif
