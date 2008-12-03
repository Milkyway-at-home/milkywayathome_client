#ifndef H_GEM_RECOMBINATION
#define H_GEM_RECOMBINATION

double* mutate(double* parent, double* min_parameters, double* max_parameters, int number_parameters);

void random_recombination(int number_parameters, double* min_parameters, double* max_parameters, double* result);

double* average_recombination(double** parents, int number_parents, int number_parameters);
double* higher_recombination(double** parents, int number_parents, int number_parameters);
double* lower_recombination(double** parents, int number_parents, int number_parameters);

double* simplex_recombination(double** parents, double* fitness, int number_parents, int number_parameters, double l1, double l2);

double* binomial_recombination(double** parents, int number_parents, int number_parameters, double crossover_rate, double crossover_scale);
double* exponential_recombination(double** parents, int number_parents, int number_parameters, double crossover_rate, double crossover_scale);

double* get_pair_sum(double **individuals, int number_individuals, int number_parameters, int number_pairs, double scale);
double* get_dir_sum(double **individuals, double *fitness, int number_individuals, int number_parameters, int number_pairs, double scale);

#endif
