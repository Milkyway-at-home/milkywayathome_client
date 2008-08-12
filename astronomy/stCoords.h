#ifndef STCOORDS_H
#define STCOORDS_H 1

void lbr2xyz( const double* lbr, double* xyz );
void xyz2lbr( const double* xyz, double* lbr );

int lbr2stream( const double* lbr, const double* spars, double* stream, int verb );
int xyz2stream( const double* xyz, const double* spars, double* stream, int verb );

void stream2lbr( const double* stream, const double* spars, double* lbr );
void stream2xyz( const double* stream, const double* spars, double* xyz );

double wedge_eta ( int wedge );
double wedge_incl( int wedge );
double get_node();

void gc2lb( int wedge, double mu, double nu, double* l, double* b );

void stripe_normal ( int wedge, double *xyz);

void lbToXyz ( double l, double b, double *xyz );

void lbToXyz(double l, double b, double *xyz); 

void xyz_mag(double* point, double offset, double* logPoint);

void xyz2lbg(double* point, double offset, double* logPoint);

void sgr_stripe_normal(int wedge, double *xyz);

#endif /* STCOORDS_H */
