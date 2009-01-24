#ifndef FGDO_REGRESSION
#define FGDO_REGRESSION

double parabolic_2d_evaluate(int number_points, double *x, double **hessian, double *gradient, double c);

void parabolic_2d_regression(int number_points, int x_length, double **x, double *y, double **hessian, double *gradient, double *c);
void parabolic_2d_regression_error(int number_points, int x_length, double **x, double *y, double **hessian, double **hessian_error, double *gradient, double *gradient_error, double *c, double *c_error);

double parabolic_center(double a, double b, double c);
void parabolic_regression(int number_points, double *x, double *y, double *a, double *b, double *c, double *a_error, double *b_error, double *c_error);

#endif
