#include <stdlib.h>

#include "settings.h"

char* working_directory = ".";

char* get_working_directory() {
	return working_directory;
}

void set_working_directory(char* wd) {
	free(working_directory);
	working_directory = wd;
}
