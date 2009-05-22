/*
 * Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
 * Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
 * and Rensselaer Polytechnic Institute.
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include <stdarg.h>
#include <stdio.h>

#include "search_log.h"
#include "../util/io_util.h"

FILE* log_open(char *search_name) {
	char logfilename[1024];
	sprintf(logfilename, "%s/%s/log", get_working_directory(), search_name);
	return fopen(logfilename, "a");
}
        
void log_fwrite_double_array(char *search_name, char* a_name, int a_length, double* a) {
	FILE *logfile = log_open(search_name);
	if (logfile == NULL) {
		printf("ERROR, COULD NOT OPEN LOGFILE");
		return;
	}
	fwrite_double_array(logfile, a_name, a_length, a);
	fclose(logfile);
}
 
void log_printf(char *search_name, char* text, ...) {
	va_list args;
	FILE *logfile = log_open(search_name);
	if (logfile == NULL) {  
		printf("ERROR, COULD NOT OPEN LOGFILE");
		return;
	}
	va_start(args, text);
	vfprintf(logfile, text, args);
	va_end(args);
	fclose(logfile);
}

FILE* error_log_open(char *search_name) {
	char logfilename[1024];
	sprintf(logfilename, "%s/%s/error", get_working_directory(), search_name);
	return fopen(logfilename, "a");
}

void error_log_printf(char *search_name, char* text, ...) {
        va_list args;
        FILE *error_file = error_log_open(search_name);
        if (error_file == NULL) {
                printf("ERROR, COULD NOT OPEN ERROR LOG FILE");
                return;
        }
        va_start(args, text);
        vfprintf(error_file, text, args);
        va_end(args);
        fclose(error_file);
}
