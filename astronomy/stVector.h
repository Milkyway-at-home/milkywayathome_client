#ifndef SPVECTOR_H
#define SPVECTOR_H 1

double** matalloc( int nrows, int ncols );
void matfree( double** mat, int nrows );

double norm( const double* vec );
void normalize( double* vec );

double dotp( const double* a, const double* b );
void crossp( const double* a, const double* b, double* prod );

double vecangle( const double* a, const double* b );

void get_transform( const double* f, const double* t, double** mat );

void do_transform( double* v, double* const* mat );

void transform_point( double* point, double** cmat, double* xsun, double* logPoint ); 

#endif /* SPVECTOR_H */
