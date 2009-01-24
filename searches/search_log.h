#ifndef FGDO_SEARCH_LOG
#define FGDO_SEARCH_LOG

#include <stdarg.h>
#include <stdio.h>

#include "search_log.h"
#include "../util/io_util.h"
#include "../util/settings.h"

FILE* log_open(char *search_name);
void log_print_double_array(char *search_name, char* a_name, int a_length, double* a);
void log_printf(char *search_name, char* text, ...);

#endif
