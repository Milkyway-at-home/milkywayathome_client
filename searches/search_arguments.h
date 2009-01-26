#ifndef FGDO_SEARCH_ARGUMENTS_H
#define FGDO_SEARCH_ARGUMENTS_H

int argument_exists(const char* name, int argc, char** argv);

int get_int_arg(const char* name, int argc, char** argv);
long get_long_arg(const char* name, int argc, char** argv);
double get_double_arg(const char* name, int argc, char** argv);

#endif
