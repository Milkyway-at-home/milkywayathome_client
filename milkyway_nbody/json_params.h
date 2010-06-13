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

#include <popt.h>
#include <json/json.h>


/* function read a named field into an enum */
/* :: (Enum a) => String -> a */
typedef int (*ReadEnum) (const char*);


/* We could use a little more than the basic json types, and maintain
 * compatability with json_type */
typedef enum
{
    nbody_type_null    = json_type_null,
    nbody_type_boolean = json_type_boolean,
    nbody_type_double  = json_type_double,
    nbody_type_int     = json_type_int,
    nbody_type_object  = json_type_object,
    nbody_type_array   = json_type_array,
    nbody_type_string  = json_type_string,

    /* new stuff */
    nbody_type_enum,   /* a string with an associated enum value */
    nbody_type_group,  /* a container */
    nbody_type_group_item /* an object with an associated type, in the container */
} nbody_type;

/* Same basic idea as a popt option table */
typedef struct _Parameter
{
    const char* name;       /* Name of the parameter */
    const nbody_type type;  /* Type of the parameter. Use json_type_int for an enum */
    const void* param;      /* Depends on type. pointer to where the value
                               should go. For an nbody_type_group, should
                               point to the type field of the group object.
                            */
    const void* dflt;       /* Default value. Use NULL to require it. Note
                               that for reading an enum from a string, use
                               the enum value and not the string name. */

    const ReadEnum conv;    /* function to read a name into an enum */
    const bool unique;      /* If this is an object, whether the parameters are supposed to be unique.
                               i.e. there's a set of options to choose from.
                            */
    const struct _Parameter* parameters; /* if an object, the subitems */
} Parameter;

/* Macros useful for making the parameter tables because they became beastly */
#define NULLPARAMETER { NULL, nbody_type_null, NULL, NULL, NULL, FALSE, NULL }

#define BASIC_PARAM(name, type, dest) { name, type, dest, NULL, NULL, FALSE, NULL }
#define DBL_PARAM(name, dest) { name, nbody_type_double, dest, NULL, NULL, FALSE, NULL }
#define INT_PARAM(name, dest) { name, nbody_type_int, dest, NULL, NULL, FALSE, NULL }
#define STR_PARAM(name, dest) { name, nbody_type_string, dest, NULL, NULL, FALSE, NULL }
#define BOOL_PARAM(name, dest) { name, nbody_type_boolean, dest, NULL, NULL, FALSE, NULL }
#define DBL_PARAM_DFLT(name, dest, dfl) { name, nbody_type_double, dest, dfl, NULL, FALSE, NULL }
#define BOOL_PARAM_DFLT(name, dest, dfl) { name, nbody_type_boolean, dest, dfl, NULL, FALSE, NULL }
#define ENUM_PARAM(name, dest, readf) { name, nbody_type_enum, dest, NULL, readf, FALSE, NULL }

/* A group where we want one of the options. dest is where the enum
 * identifier for the sub-item will go */
#define GROUP_PARAM(name, dest, items) { name, nbody_type_group, dest, NULL, NULL, TRUE, items }

/* The group item with it's enum value */
#define GROUP_PARAM_ITEM(name, val, items) { name, nbody_type_group_item, NULL, val, NULL, FALSE, items }

/* A set of parameters where all of them are required. Mostly for organization */
#define OBJ_PARAM(name, items) { name, nbody_type_object, NULL, NULL, NULL, FALSE, items }

//typedef static const Parameter PTable;


void initNBody(int argc, const char** argv);
void get_params_from_json(NBodyCtx* ctx, json_object* fileObj);
static bool warn_extra_params(json_object* obj, const char* grpName);
static void readParameterGroup(const Parameter*, json_object*, const Parameter*);

#endif /* _JSON_PARAMS_H_ */

