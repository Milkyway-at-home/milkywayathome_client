#include <math.h>
#include <string.h>
#include "atSurveyGeometry.h"
#include "stVector.h"


void atGCToEq (
	       double amu,	/* IN -- mu in degrees */
	       double anu,	/* IN -- nu in degrees */
	       double *ra,	/* OUT -- ra in degrees */
	       double *dec,	/* OUT -- dec in degrees */
	       double anode,	/* IN -- node in degrees */
	       double ainc	/* IN -- inclination in degrees */
	       )
{
  double x1, y1, z1, x2, y2, z2;
  /* Convert to radians */

//	printf("atGCToEq: amu: %g, anu: %g, anode: %g, ainc: %g\n", amu, anu, anode, ainc);

  amu = amu * at_deg2Rad;
  anu = anu * at_deg2Rad;
  anode = anode * at_deg2Rad;
  ainc = ainc * at_deg2Rad;
  /* Rotation */
  x2 = cos(amu-anode)*cos(anu);
  y2 = sin(amu-anode)*cos(anu);
  z2 = sin(anu);
  x1 = x2;
  y1 = y2*cos(ainc) - z2*sin(ainc);
  z1 = y2*sin(ainc) + z2*cos(ainc);

  *ra = atan2 (y1, x1)+anode;
  *dec = asin(z1);
  /* Convert back to degrees */
  *ra = *ra * at_rad2Deg;
  *dec = *dec * at_rad2Deg;

//  printf("atGCToEq finish: ra: %g, dec: %g\n",*ra, *dec);

  atBound2(dec, ra);
  return;
}

void atEqToGal (
		double ra,	/* IN -- ra in degrees */
		double dec,	/* IN -- dec in degrees */
		double *glong,	/* OUT -- Galactic longitude in degrees */
		double *glat	/* OUT -- Galactic latitude in degrees */
		)
{
  /* Convert to radians */
   ra = ra * at_deg2Rad;
   dec = dec * at_deg2Rad;
  /* Use SLALIB to do the actual conversion */
  slaEqgal(ra, dec, glong, glat);
  /* Convert back to degrees */
  *glong = *glong * at_rad2Deg;
  *glat = *glat * at_rad2Deg;
  atBound2(glat, glong);
  return;
}



void atBound (
	      double *angle,	/* MODIFIED -- the angle to bound in degrees*/
	      double min,	/* IN -- inclusive minimum value */
	      double max	/* IN -- exclusive maximum value */
	      ) {
  while (*angle<min) {
    *angle += 360.0;
  }
  while (*angle>=max) {
    *angle -= 360.0;
  }
  return;
}

void atBound2(
	      double *theta,	/* MODIFIED -- the -90 to 90 angle */
	      double *phi	/* MODIFIED -- the 0 to 360 angle */
	      ) {
  atBound(theta, -180.0, 180.0);
  if (fabs(*theta) > 90.0) {
    *theta = 180.0 - *theta;
    *phi += 180;
  }
  atBound(theta, -180.0, 180.0);
  atBound(phi, 0.0, 360.0);
  if (fabs(*theta)==90.0) *phi=0.;
  return;
}

void slaDcc2s ( double v[3], double *a, double *b )
{
   double x, y, z, r;
                                                                                
   x = v[0];
   y = v[1];
   z = v[2];
   r = sqrt ( x * x + y * y );
                                                                                
   *a = ( r != 0.0 ) ? atan2 ( y, x ) : 0.0;
   *b = ( z != 0.0 ) ? atan2 ( z, r ) : 0.0;
}


void slaDcs2c ( double a, double b, double v[3] )
{
   double cosb;
                                                                                
   cosb = cos ( b );
   v[0] = cos ( a ) * cosb;
   v[1] = sin ( a ) * cosb;
   v[2] = sin ( b );
}


void slaDimxv ( double dm[3][3], double va[3], double vb[3] )
{
  long i, j;
  double w, vw[3];
                                                                                
/* Inverse of matrix dm * vector va -> vector vw */
   for ( j = 0; j < 3; j++ ) {
      w = 0.0;
      for ( i = 0; i < 3; i++ ) {
         w += dm[i][j] * va[i];
      }
      vw[j] = w;
   }
                                                                                
/* Vector vw -> vector vb */
   for ( j = 0; j < 3; j++ ) {
     vb[j] = vw[j];
   }
}


void slaDmxv ( double dm[3][3], double va[3], double vb[3] )
{
   int i, j;
   double w, vw[3];
                                                                                
/* Matrix dm * vector va -> vector vw */
   for ( j = 0; j < 3; j++ ) {
      w = 0.0;
      for ( i = 0; i < 3; i++ ) {
         w += dm[j][i] * va[i];
      }
      vw[j] = w;
   }
                                                                                
/* Vector vw -> vector vb */
   for ( j = 0; j < 3; j++ ) {
      vb[j] = vw[j];
   }
}


double slaDrange ( double angle )
{
  double w;
                                                                                
  w = dmod ( angle, D2PI );
  return ( fabs ( w ) < DPI ) ? w : w - dsign ( D2PI, angle );
}


double slaDranrm ( double angle )
{
   double w;
                                                                                
   w = dmod ( angle, D2PI );
   return ( w >= 0.0 ) ? w : w + D2PI;
}


void slaEqgal ( double dr, double dd, double *dl, double *db )
{
   double v1[3], v2[3];
                                                                                
   static double rmat[3][3];
                                                                                
   rmat[0][0] = -0.054875539726;
   rmat[0][1] = -0.873437108010;
   rmat[0][2] = -0.483834985808;
   rmat[1][0] =  0.494109453312;
   rmat[1][1] = -0.444829589425;
   rmat[1][2] =  0.746982251810;
   rmat[2][0] = -0.867666135858;
   rmat[2][1] = -0.198076386122;
   rmat[2][2] =  0.455983795705;
                                                                                
/* Spherical to Cartesian */
   slaDcs2c ( dr, dd, v1 );
                                                                                
/* Equatorial to Galactic */
   slaDmxv ( rmat, v1, v2 );
                                                                                
/* Cartesian to spherical */
   slaDcc2s ( v2, dl, db );
                                                                                
/* Express in conventional ranges */
   *dl = slaDranrm ( *dl );
   *db = slaDrange ( *db );
}

