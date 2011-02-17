/* Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _JSON_PARAMS_H_
#define _JSON_PARAMS_H_

#ifdef __cplusplus
extern "C"
#endif

#include <json/json.h>
#include "nbody_types.h"

/* function read a named field into an enum */
/* :: (Enum a) => String -> a */
typedef int (*MWReadEnum) (const char*);


/* We could use a little more than the basic json types, and maintain
 * compatability with json_type */
typedef enum
{
    mw_type_null    = json_type_null,
    mw_type_boolean = json_type_boolean,
    mw_type_double  = json_type_double,
    mw_type_int     = json_type_int,
    mw_type_object  = json_type_object,
    mw_type_array   = json_type_array,
    mw_type_string  = json_type_string,

    /* new stuff */
    mw_type_vector,     /* an array of length 3 */
    mw_type_enum,       /* a string with an associated enum value */
    mw_type_one_or_many /* Try a single object, or try a list of them */
} mw_type;

/* ItemToBeRead -> NameOfItem -> JsonObject -> True on failure */
typedef mwbool (*MWReadFunc) (void*, const char*, json_object*);

/* Read a json_object into some type */
typedef mwbool (*MWArrayRead) (void*, json_object*, const void*);


/* Same basic idea as a popt option table */
typedef struct _MWParameter
{
    const char* name;       /* Name of the parameter */
    const mw_type type;     /* Type of the parameter. Use json_type_int for an enum */
    void* param;            /* Depends on type. pointer to where the value
                               should go. */
    const size_t size;      /* Size of type for a list. Should be 0 otherwise */
    unsigned int* length;   /* If an array, number of values out */
    const void* dflt;       /* Default value. Use NULL to require it. Note
                               that for reading an enum from a string, use
                               the enum value and not the string name. */

    const MWReadFunc conv; /* function to read a name into an enum. Should return -1 for failure. */
} MWParameter;

/* Terrible, terrible macros useful for making the parameter tables */
#define NULL_MWPARAMETER { NULL, mw_type_null, NULL, 0, NULL, NULL, NULL }


/* Sets of different parameters with an associated enum. e.g. a disk
 * type can have different parameters depending on what kind it is. */
typedef struct
{
    const char* name;              /* Name of the option */
    generic_enum_t enumVal;        /* Enum value correcsponding to the name */
    const MWParameter* parameters; /* Associated set of parameters */
} MWParameterSet;

#define NULL_MWPARAMETERSET { NULL, -1, NULL }


#define BASIC_PARAM(name, type, dest) { name, type, dest, 0, NULL, NULL, NULL }
#define STR_PARAM(name, dest) { name, mw_type_string, dest, 0, NULL, NULL, NULL }

#define DBL_PARAM(name, dest) { name, mw_type_double, dest, 0, NULL, NULL, NULL }
#define DBL_PARAM_DFLT(name, dest, dfl) { name, mw_type_double, dest, 0, NULL, dfl, NULL }

#define INT_PARAM(name, dest) { name, mw_type_int, dest, 0, NULL, NULL, NULL }
#define INT_PARAM_DFLT(name, dest, dfl) { name, mw_type_int, dest, 0, NULL, dfl, NULL }

#define ENUM_PARAM(name, dest, readf) { name, mw_type_enum, dest, 0, NULL, NULL, readf }
#define ENUM_PARAM_DFLT(name, dest, dflt, readf) { name, mw_type_enum, dest, 0, NULL, dflt, readf }

#define BOOL_PARAM(name, dest) { name, mw_type_boolean, dest, 0, NULL, NULL, NULL }
#define BOOL_PARAM_DFLT(name, dest, dfl) { name, mw_type_boolean, dest, 0, NULL, dfl, NULL }

#define VEC_PARAM(name, dest) { name, mw_type_vector, dest, 0, NULL, NULL, NULL }
#define VEC_PARAM_DFLT(name, dest, dfl) { name, mw_type_vector, dest, 0, NULL, dfl, NULL }

#define OBJ_PARAM(name, dest, readf) { name, mw_type_object, dest, 0, 0, NULL, readf }
#define ARRAY_PARAM(name, dest, size, length, readf) { name, mw_type_array, dest, size, length, NULL, readf }
#define ONE_OR_MANY_PARAM(name, dest, size, length, readf) { name, mw_type_one_or_many, dest, size, length, NULL, readf }


const MWParameterSet* mwReadParameterSet(const MWParameterSet* ps, const char* str, const char* name);
mwbool mwReadParameterGroup(const MWParameter* g, json_object* hdr, const char* parentName);
mwbool mwReadTypedGroup(const MWParameterSet* ps, json_object* obj, const char* name, generic_enum_t* typeRead);


#ifdef __cplusplus
}
#endif

#endif /* _JSON_PARAMS_H_ */

