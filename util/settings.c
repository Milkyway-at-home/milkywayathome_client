#include <stdlib.h>
#include <stdio.h>

#include "settings.h"

char* working_directory = ".";

char* get_working_directory() {
	return working_directory;
}

void set_working_directory(char* wd) {
	free(working_directory);
	working_directory = wd;
}

void remove_arg(int *target, int *count, char ***values) {
	int i, current;
	char** new_values;

	new_values = (char**)malloc(sizeof(char*) * ((*count) - 1));
	current = 0;
	for (i = 0; i < (*count); i++) {
		if (i == (*target)) continue;
		new_values[current] = (*values)[i];
		current++;
	}

	printf("freeing %d of %d\n", (*target), (*count));
	printf("\tvalue: %s\n", (*values)[(*target)]);
	free( (*values)[(*target)] );
	free( (*values) );

	(*values) = new_values;
	(*target)--;
	(*count)--;


}
