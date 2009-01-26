#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "search_arguments.h"


int argument_exists(const char* name, int argc, char** argv) {
	int i;
	for (i = 0; i < argc; i++) if (!strcmp(argv[i], name)) return 1;
	return 0;
}

int get_int_arg(const char* name, int argc, char** argv) {
	int i;
	for (i = 0; i < argc; i++) if (!strcmp(argv[i], name)) return atoi(argv[++i]);
	return -1;
}

long get_long_arg(const char* name, int argc, char** argv) {
	int i;
	for (i = 0; i < argc; i++) if (!strcmp(argv[i], name)) return atol(argv[++i]);
	return -1;
}

double get_double_arg(const char* name, int argc, char** argv) {
	int i;
	for (i = 0; i < argc; i++) if (!strcmp(argv[i], name)) return atof(argv[++i]);
	return -1.0;
}
