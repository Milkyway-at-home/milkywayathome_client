#ifndef TEST_ENV_UTIL_H
#define TEST_ENV_UTIL_H

#include "nbody_types.h"
#include "nbody_dwarf_potential.h"

/* Function declarations */
char* find_lua_file(const char* filename, ...);
char* find_milkyway_nbody();
int run_nbody(const char** dwarf_params, const char* lua_file);
int read_lua_parameters(const char* input_lua_file, const char** dwarf_params, 
                        real* nbody, real* nbody_baryon, 
                        Dwarf** comp1, Dwarf** comp2, real* timestep);

#endif /* TEST_ENV_UTIL_H */