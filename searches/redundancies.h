#ifndef FGDO_REDUNDANCIES_H
#define FGDO_REDUNDANCIES_H

#include "stdio.h"
#include "../util/settings.h"

#define PARAMETER_THRESHOLD 10e-10
#define FITNESS_THRESHOLD 10e-10

typedef struct redundancy {
	double fitness;
	int hostid;

	struct redundancy *next;
} REDUNDANCY; 

typedef struct redundancy_list {
	int number_redundancies;
	int number_parameters;
	double *parameters;
	char *metadata;

	REDUNDANCY *first;
	struct redundancy_list *next;
} REDUNDANCY_LIST;

typedef struct redundancies {
	int search_size;
	int number_redundancies;
	double redundancy_rate;

	REDUNDANCY_LIST *redundancy_list;
	REDUNDANCY_LIST *last_list;
} REDUNDANCIES;



void initialize_redundancies(int search_size, double redundancy_rate, REDUNDANCIES **redundancies);

int fwrite_redundancies(FILE *file, REDUNDANCIES *redundancies);
int write_redundancies(char *filename, REDUNDANCIES *redundancies);

int fread_redundancies(FILE *file, REDUNDANCIES **redundancies);
int read_redundancies(char *filename, REDUNDANCIES **redundancies);

int generate_redundancy(REDUNDANCIES *redundancies, int number_parameters, double *parameters, char *metadata);


#define REDUNDANCY_MATCH 0
#define REDUNDANCY_MISMATCH 1
#define REDUNDANCY_DUPLICATE 2
#define REDUNDANCY_NOT_EQUAL 3
int redundancy_match(REDUNDANCY_LIST *current, int number_parameters, double fitness, double *parameters, int hostid);

#define VERIFY_VALID 0
#define VERIFY_IN_PROGRESS 1
#define VERIFY_INSERT 2
#define VERIFY_NOT_INSERTED 3
int verify_with_insert(REDUNDANCIES *redundancies, int number_parameters, double fitness, double *parameters, char *metadata, int hostid);
int verify_without_insert(REDUNDANCIES *redundancies, int number_parameters, double fitness, double *parameters, char *metadata, int hostid);

#endif
