#include "math.h"
#include "stdlib.h"
#include "stdio.h"

#include "redundancy.h"

#define FITNESS_MATCH_TOLERANCE 10e-10
#define PARAMETERS_MATCH_TOLERANCE 10e-15

int fitness_match(double f1, double f2) {
        return fabs(f1 - f2) < FITNESS_MATCH_TOLERANCE;
}

int parameters_match(int n, double *p1, double *p2) {
	int i;
        for (i = 0; i < n; i++) {
                if (fabs(p1[i] - p2[i]) > PARAMETERS_MATCH_TOLERANCE) return 0;
        }
        return 1;
}

void free_redundancy(REDUNDANCY **redundancy) {
	REDUNDANCY *r = (*redundancy);
	free(r->parameters);
	free(r->velocity);
	r->next = NULL;
	free(r);
	(*redundancy) = NULL;
}

void new_redundancy(REDUNDANCY **redundancy, double fitness, int n, double *parameters, double *velocity) {
	int i;
	REDUNDANCY *r;
	(*redundancy) = (REDUNDANCY*)malloc(sizeof(REDUNDANCY));
	r = (*redundancy);
	r->fitness = fitness;
	r->parameters = (double*)malloc(sizeof(double) * n);
	r->velocity = (double*)malloc(sizeof(double) * n);
	for (i = 0; i < n; i++) {
		r->parameters[i] = parameters[i];
		r->velocity[i] = velocity[i];
	}
	r->next = NULL;
}

int fread_redundancy(FILE *file, REDUNDANCY **r, int n, int *particle) {
	double fitness, *parameters, *velocity;
	int i;
	if (2 != fscanf(file, "[%d] fitness: %lf, parameters:", particle, &fitness)) return 0;
	parameters = (double*)malloc(sizeof(double) * n);
	for (i = 0; i < n; i++) {
		if (1 != fscanf(file, " %lf", &(parameters[i]))) {
			free(parameters);
			(*r) = NULL;
			return 0;
		}
	}
	fscanf(file, ", velocity:");
	velocity = (double*)malloc(sizeof(double) * n);
	for (i = 0; i < n; i++) {
		if (1 != fscanf(file, " %lf", &(velocity[i]))) {
			free(parameters);
			free(velocity);
			(*r) = NULL;
			return 0;
		}
	}
	fscanf(file, "\n");

	(*r) = (REDUNDANCY*)malloc(sizeof(REDUNDANCY));
	(*r)->fitness = fitness;
	(*r)->velocity = velocity;
	(*r)->parameters = parameters;
	(*r)->next = NULL;
	return 1;
}


int fwrite_redundancy(FILE *file, REDUNDANCY *r, int n, int particle) {
	int i;
	fprintf(file, "[%d] fitness: %.20lf, parameters:", particle, r->fitness);
	for (i = 0; i < n; i++) fprintf(file, " %.20lf", r->parameters[i]);
	fprintf(file, ", velocity:");
	for (i = 0; i < n; i++) fprintf(file, " %.20lf", r->velocity[i]);
	fprintf(file, "\n");
	return 1;
}

int check_redundancy(REDUNDANCY *search_current, REDUNDANCY *r, int point, SEARCH_PARAMETERS *sp) {
        /**
	 ** See if this particle has any saved local best values
	 ** If this matches a local best value update the local best
	 **      remove the match from the queue
	 **      remove all queued matches with lower fitness
	 ** If this matches parameters but not fitness, remove the match from the queue
	 **/
        REDUNDANCY *r, *r_prev, *n, *n_prev;
        int match;
        r = sx->redundancies[point];
        r_prev = NULL;
        while (r != NULL) {
                if (parameters_match(sx->number_parameters, r->parameters, sp->parameters)) {
                        match = fitness_match(r->fitness, sp->fitness);

                        if (match) {
                                n_prev = NULL;
                                n = sx->redundancies[point];
				while (n != NULL) {
					if (n->fitness <= sp->fitness) {
						REDUNDANCY *temp;
						temp = n->next;
						if (n_prev == NULL) {
							sx->redundancies[point] = n->next;
						} else {
							n_prev->next = n->next;
						}
						free_redundancy(&n);
						n = temp;
					} else {
						n_prev = n;
						n = n->next;
					}
				}
			} else {
				if (r_prev == NULL) {
					sx->redundancies[point] = r->next;
					free_redundancy(&r);
				} else {
					r_prev->next = r->next;
					free_redundancy(&r);
				}
			}
			sx->current_redundancy = NULL;
			return match;
		}
		r_prev = r;
		r = r->next;
	}
	return 0;
}
