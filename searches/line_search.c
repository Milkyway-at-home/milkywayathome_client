#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "../evaluation/evaluator.h"
#include "../util/io_util.h"

/********
	*	Stopping Conditions
 ********/
int		nquad = 3;
double		tol = 1e-6;

double evaluate_step(double* point, double step, double* direction, int number_parameters, double* current_point) {
	int i;
	for (i = 0; i < number_parameters; i++) {
		current_point[i] = point[i] - (step * direction[i]);
	}
	return evaluate(current_point);
}

int synchronous_line_search(double* point, double initial_fitness, double* direction, int number_parameters, double** minimum, double* fitness) {
	int i, jump, evaluations_done;
	double f1, f2, f3, fs;
	double d1, d2, d3, dstar;
	double step, a, b;
	double* current_point;

	current_point = (double*)malloc(sizeof(double) * number_parameters);

	printf("\tline search starting at fitness: %lf\n", initial_fitness);
	print_double_array(stdout, "\tinitial point: ", number_parameters, point);
	/********
		*	Find f1 > f2
	 ********/
	step = 1.0;
	f1 = initial_fitness;
	// f2 = evaluate( point + (direction * step) );
	f2 = evaluate_step(point, step, direction, number_parameters, current_point);
	evaluations_done = 1;
	printf("\t\tloop 1, evaluations: %d, step: %lf, fitness: %lf\n", evaluations_done, step, f2);

	while (!isnan(f2) && f1 < f2 && step > tol) {
		step /= 2.0;

		// f2 = evaluate( point + (direction * step) );
		f2 = evaluate_step(point, step, direction, number_parameters, current_point);
		evaluations_done++;
		printf("\t\tloop 1, evaluations: %d, step: %lf, fitness: %lf\n", evaluations_done, step, f2);
	}
	if (isnan(f2)) return -1;
	if (step < tol) return -1;

	/********
		*	Find f1 > f2 < f3
	 ********/
	jump = 2.0;
	d1 = 0.0;
	d2 = 1.0;
	d3 = 2.0;
	// f3 = evaluate( point = (d3 * step * direction) );
	f3 = evaluate_step(point, d3 * step, direction, number_parameters, current_point);
	evaluations_done++;
	printf("\t\tloop 2, evaluations: %d, step: %lf, fitness: %lf\n", evaluations_done, d3 * step, f3);

	while (f3 <= f2 && !isnan(f1 - f2 - f3)) {
		d1 = d2;
		f1 = f2;
		d2 = d3;
		f2 = f3;
		d3 = jump * d3;

		// f3 = evaluate( point + (d3 * step * direction) );
		f3 = evaluate_step(point, d3 * step, direction, number_parameters, current_point);
		evaluations_done++;
		printf("\t\tloop 2, evaluations: %d, step: %lf, fitness: %lf\n", evaluations_done, d3 * step, f3);
	}
	if (isnan(f1 - f2 - f3)) return -1;

	/********
		*	Find minimum
	 ********/
	fs = 0.0;
	for (i = 0; i < nquad; i++) {
		a = (d1*d1)*(f2-f3) + (d2*d2)*(f3-f1) + (d3*d3)*(f1-f2);
		b = d1*(f2-f3) + d2*(f3-f1) + d3*(f1-f2);
		dstar = 0.5 * (a/b);

		if (dstar < d2 + tol && dstar > d2 - tol) dstar = d2 + tol;

		// fs = evaluate(point + (dstar * step * direction));
		fs = evaluate_step(point, dstar * step, direction, number_parameters, current_point);
		evaluations_done++;
		printf("\t\tloop 3, evaluations: %d, step: %lf, fitness: %lf\n", evaluations_done, (dstar * step), fs);
		if (isnan(fs)) return -1;

		if (dstar > d2 + tol) {
			if (fs > f2) {
				d3 = dstar;
				f3 = fs;
			} else {
				d1 = d2;
				f1 = f2;
				d2 = dstar;
				f2 = fs;
			}
		} else {
			if (fs > f2) {
				d1 = dstar;
				f1 = fs;
			} else {
				d3 = d2;
				f3 = f2;
				d2 = dstar;
				f2 = fs;
			}
		}
	}

	(*minimum) = (double*)malloc(sizeof(double) * number_parameters); 
	memcpy((*minimum), current_point, sizeof(double) * number_parameters);
	(*fitness) = fs;

	free(current_point);
	return evaluations_done;
}
