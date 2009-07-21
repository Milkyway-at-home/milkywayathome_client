#include "redundancies.h"
#include "../util/settings.h"

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

void new_redundancy(double fitness, int hostid, REDUNDANCY **redundancy) {
	(*redundancy) = (REDUNDANCY*)malloc(sizeof(REDUNDANCY));
	(*redundancy)->fitness = fitness;
	(*redundancy)->hostid = hostid;
	(*redundancy)->next = NULL;
}

void free_redundancy_list(REDUNDANCY_LIST **redundancy_list) {
	REDUNDANCY *current_redundancy, *temp;

	current_redundancy = (*redundancy_list)->first;
	while (current_redundancy != NULL) {
		temp = current_redundancy;
		current_redundancy = current_redundancy->next;
		free(temp);
	}

	free( (*redundancy_list)->parameters );
	free( (*redundancy_list)->metadata );
	free( (*redundancy_list) );
}

void new_redundancy_list(int number_parameters, double fitness, double *parameters, char *metadata, int hostid, REDUNDANCY_LIST **redundancy_list) {
	int i;
	(*redundancy_list) = (REDUNDANCY_LIST*)malloc(sizeof(REDUNDANCY_LIST));

	(*redundancy_list)->number_parameters = number_parameters;
	(*redundancy_list)->parameters = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) (*redundancy_list)->parameters[i] = parameters[i];

	(*redundancy_list)->metadata = (char*)malloc(sizeof(char) * METADATA_SIZE);
	memcpy((*redundancy_list)->metadata, metadata, sizeof(char) * METADATA_SIZE);

	(*redundancy_list)->next = NULL;

	new_redundancy(fitness, hostid, &((*redundancy_list)->first));
}

void initialize_redundancies(REDUNDANCIES **redundancies) {
	(*redundancies) = (REDUNDANCIES*)malloc(sizeof(REDUNDANCIES));
	(*redundancies)->redundancy_list = NULL;
	(*redundancies)->last_list = NULL;
}


int fwrite_redundancies(FILE *file, REDUNDANCIES *redundancies) {
	int i;
	REDUNDANCY_LIST *current_list;
	REDUNDANCY *current_redundancy;

	current_list = redundancies->redundancy_list;
	fprintf(file, "redundancies:\n");
	while (current_list != NULL) {
		fprintf(file, "[%d] [%.20lf", current_list->number_parameters, current_list->parameters[0]);
		for (i = 1; i < current_list->number_parameters; i++) fprintf(file, " %.20lf", current_list->parameters[i]);
		fprintf(file, "] [%s] ", current_list->metadata);

		current_redundancy = current_list->first;
		while (current_redundancy != NULL) {
			fprintf(file, "[%.20lf %d]", current_redundancy->fitness, current_redundancy->hostid);
			current_redundancy = current_redundancy->next;
		}
		fprintf(file, "\n");
		current_list = current_list->next;
	}
	fprintf(file, "finished\n");
	return 0;
}

int fread_redundancy(FILE *file, REDUNDANCY **redundancy) {
	(*redundancy) = (REDUNDANCY*)malloc(sizeof(REDUNDANCY));

	if (2 == fscanf(file, "[%lf %d]", &((*redundancy)->fitness), &((*redundancy)->hostid))) {
		return 1;
	} else {
		free(*redundancy);
		*redundancy = NULL;
		return 0;
	}
}

int fread_redundancy_list(FILE *file, REDUNDANCY_LIST **redundancy_list) {
	int i;
	char c;
	REDUNDANCY *current_redundancy;

	(*redundancy_list) = (REDUNDANCY_LIST*)malloc(sizeof(REDUNDANCY_LIST));

	if (1 != fscanf(file, "[%d] [", &((*redundancy_list)->number_parameters))) {
		free_redundancy_list(redundancy_list);
		(*redundancy_list) = NULL;
		return 0;
	}

	(*redundancy_list)->parameters = (double*)malloc(sizeof(double) * (*redundancy_list)->number_parameters);
	for (i = 0; i < (*redundancy_list)->number_parameters; i++) {
		fscanf(file, "%lf", &((*redundancy_list)->parameters[i]));
		fgetc(file);
	}
	fscanf(file, " [");

	(*redundancy_list)->metadata = (char*)malloc(sizeof(char) * METADATA_SIZE);
	c = fgetc(file);
	while (c != ']' && i < METADATA_SIZE) {
		(*redundancy_list)->metadata[i] = c;
		i++;
		c = fgetc(file);
	}
	fgetc(file);

	fread_redundancy(file, &current_redundancy);
	(*redundancy_list)->first = current_redundancy;
	while (current_redundancy != NULL) {
		fread_redundancy(file, &(current_redundancy->next));
		current_redundancy = current_redundancy->next;
	}
	fscanf(file, "\n");
	return 0;
}


int fread_redundancies(FILE *file, REDUNDANCIES **redundancies) {
	REDUNDANCY_LIST *current_list;

	initialize_redundancies(redundancies);

	fscanf(file, "redundancies:\n");

	fread_redundancy_list(file, &current_list);
	(*redundancies)->redundancy_list = current_list;

	while (current_list != NULL) {
		fread_redundancy_list(file, &(current_list->next));
		current_list = current_list->next;
	}

	fscanf(file, "finished\n");

	return 0;
}

int read_redundancies(char *filename, REDUNDANCIES **redundancies) {
	FILE *file = fopen(filename, "r");
	if (file == NULL) return -1;

	fread_redundancies(file, redundancies);
	fclose(file);
	return 0;
}

int write_redundancies(char *filename, REDUNDANCIES *redundancies) {
	FILE *file = fopen(filename, "w");
	if (file == NULL) return -1;

	fwrite_redundancies(file, redundancies);
	fclose(file);
	return 0;
}


void generate_redundancy(REDUNDANCIES *redundancies, int number_parameters, double *parameters, char *metadata) {
	int i;
	if (redundancies->redundancy_list == NULL) return;

	for (i = 0; i < number_parameters; i++) parameters[i] = redundancies->redundancy_list->parameters[i];
	memcpy(metadata, redundancies->redundancy_list->metadata, sizeof(char) * METADATA_SIZE);

	redundancies->last_list->next = redundancies->redundancy_list;
	redundancies->redundancy_list = redundancies->redundancy_list->next;
	redundancies->last_list = redundancies->last_list->next;
	redundancies->last_list->next = NULL;
}

void append_redundancy(REDUNDANCY_LIST *redundancy_list, double fitness, int hostid) {
	REDUNDANCY *current;

	if (redundancy_list->first == NULL) {
		new_redundancy(fitness, hostid, &(redundancy_list->first));
		return;
	}
	current = redundancy_list->first;
	while (current->next != NULL) current = current->next;
	new_redundancy(fitness, hostid, &(current->next));
}


int redundancy_match(REDUNDANCY_LIST *current, int number_parameters, double fitness, double *parameters, int hostid) {
	int i;
	REDUNDANCY *current_redundancy;

	if (number_parameters != current->number_parameters) return REDUNDANCY_NOT_EQUAL;
	
	for (i = 0; i < current->number_parameters; i++) {
//		printf("comparing parameters[%d]: %.20lf %.20lf\n", i, parameters[i], current->parameters[i]);
		if (fabs(parameters[i] - current->parameters[i]) > PARAMETER_THRESHOLD) {
			return REDUNDANCY_NOT_EQUAL;
		}
	}

	current_redundancy = current->first;
	while (current_redundancy != NULL) {
//		printf("comparing hostid: %d %d, fitness: %.20lf %.20lf\n", hostid, current_redundancy->hostid, fitness, current_redundancy->fitness);
		if (hostid == current_redundancy->hostid) return REDUNDANCY_DUPLICATE;
		else if (fabs(fitness - current_redundancy->fitness) < FITNESS_THRESHOLD) {
			return REDUNDANCY_MATCH;
		}
	}

	return REDUNDANCY_MISMATCH;
}


int verify_with_insert(REDUNDANCIES *redundancies, int number_parameters, double fitness, double *parameters, char *metadata, int hostid) {
	REDUNDANCY_LIST *current;
	REDUNDANCY_LIST *previous;
	REDUNDANCY_LIST *new_redundancy;

	current = redundancies->redundancy_list;
	if (current == NULL) {
		new_redundancy_list(number_parameters, fitness, parameters, metadata, hostid, &(redundancies->redundancy_list));
		redundancies->last_list = redundancies->redundancy_list;
		return VERIFY_INSERT;
	}

	previous = NULL;

	while (current != NULL) {
		int match = redundancy_match(current, number_parameters, fitness, parameters, hostid);
		if (match == REDUNDANCY_MATCH) {
			/**
			 * Remove the redundancy
			 */
			if (previous == NULL) redundancies->redundancy_list = current->next;
			else previous->next = current->next;

			if (current->next == NULL) redundancies->last_list = previous;
			free_redundancy_list(&current);

			return VERIFY_VALID;
		} else if (match == REDUNDANCY_MISMATCH) {
			/**
			 * Append to the redundancy list
			 */
			append_redundancy(current, fitness, hostid);

			return VERIFY_IN_PROGRESS;
		} else if (match == REDUNDANCY_DUPLICATE) {
			/**
			 * do nothing 
			 */
			return VERIFY_IN_PROGRESS;
		} else if (match == REDUNDANCY_NOT_EQUAL) {
			/**
			 *	Keep iterating
			 */
		}
		previous = current;
		current = current->next;
	}

	new_redundancy_list(number_parameters, fitness, parameters, metadata, hostid, &(new_redundancy));
	previous->next = new_redundancy;
	redundancies->last_list = new_redundancy;

	return VERIFY_INSERT;
}

int verify_without_insert(REDUNDANCIES *redundancies, int number_parameters, double fitness, double *parameters, char *metadata, int hostid) {
	REDUNDANCY_LIST *current;
	REDUNDANCY_LIST *previous;
	
	previous = NULL;
	current = redundancies->redundancy_list;

	while (current != NULL) {
		int match = redundancy_match(current, number_parameters, fitness, parameters, hostid);
		if (match == REDUNDANCY_MATCH) {
			/**
			 * Remove the redundancy
			 */
			if (previous == NULL) redundancies->redundancy_list = current->next;
			else previous->next = current->next;

			if (current->next == NULL) redundancies->last_list = previous;
			free_redundancy_list(&current);

			return VERIFY_VALID;
		} else if (match == REDUNDANCY_MISMATCH) {
			/**
			 * Append to the redundancy list
			 */
			append_redundancy(current, fitness, hostid);

			return VERIFY_IN_PROGRESS;
		} else if (match == REDUNDANCY_DUPLICATE) {
			/**
			 * do nothing 
			 */
			return VERIFY_IN_PROGRESS;
		} else if (match == REDUNDANCY_NOT_EQUAL) {
			/**
			 *	Keep iterating
			 */
		}
		previous = current;
		current = current->next;
	}

	return VERIFY_NOT_INSERTED;
}
