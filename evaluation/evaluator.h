#ifndef GEM_EVALUATOR_H
#define GEM_EVALUATOR_H

double (*evaluate)(double*);

void evaluator__init(int *number_arguments, char*** arguments, void (*read_data)(int, int));
void evaluator__init_integral(void (*i_f)(double*, double**), int i_p_l, void (*i_c)(double*, int, double**), int i_r_l);
void evaluator__init_likelihood(void (*l_f)(double*, double**), int l_p_l, double (*l_c)(double*, int), int l_r_l);

#endif
