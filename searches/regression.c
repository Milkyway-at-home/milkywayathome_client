#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "regression.h"

#include "../util/matrix.h"
#include "../util/io_util.h"

void parabolic_2d_helper(int number_points, int x_length, double **x, double *y, double **hessian, double *gradient, double *c, double ***Xni, int *X_length) {
	int i, j, k, current, pos;
	double **X;
	double **Xt;		//X transpose
	double **Xn;		//X norm
	double **XniXt;		//X norm inverted * X transpose
	double *B;		//results

	(*X_length) = 1 + x_length + x_length;
	for (i = x_length - 1; i > 0; i--) (*X_length) += i;

	X = (double**)malloc(sizeof(double*) * number_points);
	for (i = 0; i < number_points; i++) {
		X[i] = (double*)malloc(sizeof(double) * (*X_length));
		X[i][0] = 1;
		for (j = 0; j < x_length; j++) {
			X[i][j + 1] = x[i][j];
			X[i][j + x_length + 1] = 0.5 * x[i][j] * x[i][j];
		}
		current = 0;
		for (j = 0; j < x_length; j++) {
			for (k = j+1; k < x_length; k++) {
				X[i][1 + x_length + x_length + current] = x[i][j] * x[i][k];
				current++;
			}
		}
	}

	///B = (Xt * X)^-1 * Xt * Y
	matrix_transpose(X, number_points, (*X_length), &Xt);
	matrix_multiply(Xt, (*X_length), number_points, X, number_points, (*X_length), &Xn);
	matrix_invert(Xn, (*X_length), (*X_length), Xni);
	matrix_multiply((*Xni), (*X_length), (*X_length), Xt, (*X_length), number_points, &XniXt);
	matrix_vector_multiply(XniXt, (*X_length), number_points, y, number_points, &B);
	(*c) = B[0];

	for (i = 0; i < x_length; i++) {
		gradient[i] = B[1 + i];
		hessian[i][i] = B[1 + x_length + i];
	}

	current = 0;
	for (j = 0; j < x_length; j++) {
		for (k = j+1; k < x_length; k++) {
			pos = 1 + x_length + x_length + current;
			hessian[j][k] = B[pos];
			hessian[k][j] = B[pos];
			current++;
		}
	}

	/********
		*	free X, Xt, Xn, Xni, XniXt
	 ********/
	for (i = 0; i < (*X_length); i++) {
		free(Xt[i]);
		free(Xn[i]);
		free(XniXt[i]);
	}
	for (i = 0; i < number_points; i++) {
		free(X[i]);
	}
	free(X);
	free(Xn);
	free(XniXt);
	free(Xt);
	free(B);
}

void parabolic_2d_regression(int number_points, int x_length, double **x, double *y, double **hessian, double *gradient, double *c) {
	int i, X_length;
	double **Xni;		//X norm inverted
	parabolic_2d_helper(number_points, x_length, x, y, hessian, gradient, c, &Xni, &X_length);

	for (i = 0; i < X_length; i++) free(Xni[i]);
	free(Xni);
}

double parabolic_2d_evaluate(int x_length, double *x, double **hessian, double *gradient, double c) {
	int j, k;
	double temp = c;
	for (j = 0; j < x_length; j++) temp += gradient[j] * x[j];
	for (j = 0; j < x_length; j++) {
		for (k = 0; k < x_length; k++) {
			temp += hessian[k][j] * x[k] * x[k];
		}
	}
	return temp;
}

void parabolic_2d_regression_error(int number_points, int x_length, double **x, double *y, double **hessian, double **hessian_error, double *gradient, double *gradient_error, double *c, double *c_error) {
	int i, j, k, current, pos, X_length;
	double **Xni;		//X norm inverted
	double S, temp;

	parabolic_2d_helper(number_points, x_length, x, y, hessian, gradient, c, &Xni, &X_length);

	S = 0;
	for (i = 0; i < number_points; i++) {
		temp = y[i] - parabolic_2d_evaluate(x_length, x[i], hessian, gradient, (*c));
		S += temp * temp;
	}

	(*c_error) = sqrt((S * Xni[0][0]) / (number_points - X_length));
	for (i = 0; i < x_length; i++) {
		gradient_error[i] = sqrt((S * Xni[i + 1][i + 1]) / (number_points - X_length));
		hessian_error[i][i] = sqrt((S * Xni[i + x_length + 1][i + x_length + 1]) / (number_points - X_length));
	}

	current = 0;
	for (j = 0; j < x_length; j++) {
		for (k = j+1; k < x_length; k++) {
			pos = 1 + x_length + x_length + current;
			hessian_error[j][k] = sqrt((S * Xni[pos][pos]) / (number_points - X_length));
			hessian_error[k][j] = hessian_error[j][k];
			current++;
		}
	}
	for (i = 0; i < X_length; i++) free(Xni[i]);
	free(Xni);
}


double parabolic_center(double a, double b, double c) {
	return -(b/(2*a));
}

void parabolic_regression(int number_points, double *x, double *y, double *a, double *b, double *c, double *a_error, double *b_error, double *c_error) {
	int i, j, k;
	double **X, **X_norm, *rhs, *B, **X_norm_inverse;
	double S, temp;
	int p;

	rhs = (double*)malloc(sizeof(double) * 3);
	X_norm = (double**)malloc(sizeof(double*) * 3);
	for (i = 0; i < 3; i++) X_norm[i] = (double*)malloc(sizeof(double) * 3);

	X = (double**)malloc(sizeof(double*) * number_points);
	for (i = 0; i < number_points; i++) {
		X[i] = (double*)malloc(sizeof(double) * 3);
		X[i][0] = 1;
		X[i][1] = x[i];
		X[i][2] = x[i] * x[i];
	}
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			X_norm[i][j] = 0;
			for (k = 0; k < number_points; k++) {
				X_norm[i][j] += X[k][i] * X[k][j];
			}
		}
	}
	for (i = 0; i < 3; i++) {
		rhs[i] = 0;
		for (j = 0; j < number_points; j++) {
			rhs[i] += X[j][i] * y[j];
		}
	}

	matrix_invert(X_norm, 3, 3, &X_norm_inverse);
	matrix_vector_multiply(X_norm_inverse, 3, 3, rhs, 3, &B);

	printf("B: %.20lf, %.20lf, %.20lf\n", B[2], B[1], B[0]);

	(*a) = B[2];
	(*b) = B[1];
	(*c) = B[0];

	S = 0;
	for (i = 0; i < number_points; i++) {
		temp = y[i] - (x[i] * x[i] * (*a) + x[i] * (*b) + (*c));
		S += temp * temp;
	}

	p = 3;
	(*a_error) = sqrt( (S * X_norm_inverse[2][2]) / (number_points - p) );
	(*b_error) = sqrt( (S * X_norm_inverse[1][1]) / (number_points - p) );
	(*c_error) = sqrt( (S * X_norm_inverse[0][0]) / (number_points - p) );

	for (i = 0; i < number_points; i++) free(X[i]);
	free(X);
	free(B);
	for (i = 0; i < 3; i++) free(X_norm_inverse[i]);
	free(X_norm_inverse);
}
