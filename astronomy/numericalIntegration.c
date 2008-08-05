#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define EPS 3.0e-11
#define PI 3.1415926535897932384626433832795028841971693993751

/*Gauss-Legendre quadrature taken from Numerical Recipes in C*/
void gaussLegendre(double x1, double x2, double x[], double w[], int n)
{
	int m, j, i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	m = (n+1)/2;
	xm = 0.5*(x2+x1);
	xl = 0.5*(x2-x1);
	
	//fprintf(stderr, "m = %d: xm = %g: xl = %g\n", m, xm, xl);
	for (i=1; i<=m; i++) {
		//fprintf(stderr, "starting iteration %d of outer loop\n", i);
		z = cos(PI*(i-0.25)/(n+0.5));
		do {
			p1 = 1.0;
			p2 = 0.0;
			for (j = 1; j <= n; j++) {
				//fprintf(stderr, "starting iteration %d of inner loop\n", i);
				p3 = p2;
				p2 = p1;
				p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp = n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z = z1 - p1/pp;
			//fprintf(stderr, "z-z1 = %g\n", fabs(z-z1));
		} while (fabs(z-z1) > EPS);
		
		x[i-1] = xm - xl *z;
		x[n-i] = xm + xl *z;
		w[i-1] = 2.0*xl/((1.0-z*z)*pp*pp);
		w[n-i] = w[i-1];
	}
}

double qgaus(double (*func)(double, int), double a, double b, int n, int wedge) {
	int j;
	double xr, xm, dx, s;
	double *x;
	double *w;

	x = (double*)malloc(sizeof(double)*n);
	w = (double*)malloc(sizeof(double)*n);

	gaussLegendre(-1.0, 1.0, x, w, n);

	xm = 0.5*(b+a);
	xr = 0.5*(b-a);
	s = 0;
	for (j = 0; j < n; j++) {
		dx = xr*x[j];
		s += w[j]*((*func)(xm+dx, wedge));
	}

	free(x);
	free(w);
	return s *= xr;
}

double qgaus_stream(double (*func)(double, int, int), double a, double b, int n, int wedge, int sgr_coordinates) {
        int j;
        double xr, xm, dx, s;
        double *x;
        double *w;

        x = (double*)malloc(sizeof(double)*n);
        w = (double*)malloc(sizeof(double)*n);

        gaussLegendre(-1.0, 1.0, x, w, n);

        xm = 0.5*(b+a);
        xr = 0.5*(b-a);
        s = 0;
        for (j = 0; j < n; j++) {
                dx = xr*x[j];
                s += w[j]*((*func)(xm+dx, wedge, sgr_coordinates));
        }

        free(x);
        free(w);
        return s *= xr;
}
