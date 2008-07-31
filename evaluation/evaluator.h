#ifndef GEM_EVALUATOR_H
#define GEM_EVALUATOR_H

double (*evaluate)(double*);

void evaluator__init_data(void (*read_data)(int, int));
void evaluator__init_integral(double* (*i_f)(double*), int i_p_l, double* (*i_c)(double*, int), int i_r_l);
void evaluator__init_likelihood(double* (*l_f)(double*), int l_p_l, double (*l_c)(double*, int), int l_r_l);

#endif
