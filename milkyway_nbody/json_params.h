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
typedef int (*StrToEnum) (const char*);

/* Same basic idea as a popt option table */
typedef struct
{
    const char* name; /* Name of the parameter */
    json_type type;   /* Type of the parameter. Use json_type_int for an enum */
    void* param;      /* Depends on type. pointer to where the value should go */
    StrToEnum conv;   /* function to read a name into an enum */

    void* dflt;       /* Default value. Use NULL to require it. Note
                         that for reading an enum from a string, use
                         the enum value and not the string name. */
} Parameter;

typedef struct
{
    const char* name;
    const Parameter* parameters;
} ParameterGroup;


void initNBody(int argc, const char** argv);
void get_params_from_json(json_object* fileObj);

#endif /* _JSON_PARAMS_H_ */

