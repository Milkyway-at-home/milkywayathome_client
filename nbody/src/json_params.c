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

#include <string.h>
#include <assert.h>
#include "nbody_config.h"
#include "json_params.h"
#include "nbody_priv.h"
#include "milkyway_util.h"
#include "milkyway_extra.h"

static nbody_type json_object_get_type_safe(json_object* obj)
{
    return obj ? json_object_get_type(obj) : nbody_type_null;
}

static bool json_object_is_number(json_object* obj)
{
    return obj && (json_object_is_type(obj, json_type_double) || json_object_is_type(obj, json_type_int));
}

/* also works for json_object_type */
static const char* showNBodyType(nbody_type bt)
{
    switch (bt)
    {
        case nbody_type_null:
            return "null";
        case nbody_type_boolean:
            return "bool";
        case nbody_type_double:
            return "double";
        case nbody_type_int:
            return "int";
        case nbody_type_object:
            return "object";
        case nbody_type_array:
            return "array";
        case nbody_type_vector:
            return "vector";
        case nbody_type_string:
            return "string";
        case nbody_type_enum:
            return "enum";
        default:
            warn("Trying to show unknown nbody_type %d\n", bt);
            return "<unknown type>";
    }
}

static void typeWarning(nbody_type expected, const char* name, const char* pname, json_object* obj)
{
    warn("Error: expected type %s for '%s' in '%s', but got %s\n",
         showNBodyType(expected),
         name,
         pname,
         showNBodyType(json_object_get_type_safe(obj)));
}

static bool checkIsObject(json_object* obj, const char* name)
{
    nbody_type readType;

    readType = json_object_get_type(obj);
    if (readType != nbody_type_object)
    {
        warn("Got wrong type '%s' for expected object '%s'\n", showNBodyType(readType), name);
        return TRUE;
    }

    return FALSE;
}

/* Iterate through remaining keys in obj to provide useful warnings of
 * unknown parameters in the file. Returns true if any found. */
static bool warnExtraParams(json_object* obj, const char* grpName)
{
    bool haveExtra = FALSE;

    json_object_object_foreach(obj,key,val)
    {
        haveExtra = TRUE;
        warn("Warning: In group '%s': Unknown field '%s': '%s'\n",
             grpName,
             key,
             json_object_to_json_string(val));
    }

    return haveExtra;
}

/* From a given set of options with different parameters, find the matching set */
const ParameterSet* readParameterSet(const ParameterSet* ps, const char* str, const char* name)
{
    const ParameterSet* p;

    p = ps;
    while (p->name && strcasecmp(str, p->name))
        ++p;

    if (!p->name)  /* Report available types on failure */
    {
        warn("Didn't get type of group '%s':\nOptions are:\n", name);
        p = ps;
        while (p->name)
        {
            warn("\t%s\n", p->name);
            ++p;
        }
        return NULL;
    }

    return p;
}

static bool readDouble(const Parameter* p, const char* pname, json_object* obj, bool useDflt)
{
    /* json_type_int and double are OK for numbers. i.e. you can leave off the decimal.
       We don't want the other conversions, which just give you 0.0 for anything else.
    */

    if (useDflt)
        *((real*) p->param) = *((real*) p->dflt);
    else if (json_object_is_number(obj))
    {
        *((real*) p->param) = (real) json_object_get_double(obj);
    }
    else
    {
        typeWarning(nbody_type_double, p->name, pname, obj);
        return TRUE;
    }

    return FALSE;
}

static bool readInt(const Parameter* p, const char* pname, json_object* obj, bool useDflt)
{
    /* I don't think any of the conversions are acceptable */
    if (useDflt)
        *((int*) p->param) = *((int*) p->dflt);
    else if (json_object_is_type(obj, json_type_int))
        *((int*) p->param) = json_object_get_int(obj);
    else
    {
        typeWarning(nbody_type_int, p->name, pname, obj);
        return TRUE;
    }

    return FALSE;
}

static bool readBool(const Parameter* p, const char* pname, json_object* obj, bool useDflt)
{
    /* CHECKME: Size */
    if (useDflt)
        *((bool*) p->param) = *((int*) p->dflt);
    else if (json_object_is_type(obj, json_type_boolean))
        *((bool*) p->param) = (bool) json_object_get_boolean(obj);
    else
    {
        typeWarning(nbody_type_boolean, p->name, pname, obj);
        return TRUE;
    }

    return FALSE;
}

static bool readString(const Parameter* p, const char* pname, json_object* obj, bool useDflt)
{
    /* The json_object has ownership of the string so we need to copy it. */
    if (useDflt)
        *((char**) p->param) = strdup(*((char**) p->dflt));
    else if (json_object_is_type(obj, json_type_string))
    {
        *((char**) p->param) = strdup(json_object_get_string(obj));
    }
    else
    {
        typeWarning(nbody_type_string, p->name, pname, obj);
        return TRUE;
    }

    return FALSE;
}

static bool readVector(const Parameter* p, const char* pname, json_object* obj, bool useDflt)
{
    int i, arrLen;
    array_list* arr;
    json_object* tmp;

    /* FIXME: Right now assuming no default vectors, groups etc. will be used */

    if (!json_object_is_type(obj, json_type_array))
    {
        typeWarning(nbody_type_vector, p->name, pname, obj);
        return TRUE;
    }

    arr = json_object_get_array(obj);
    arrLen = json_object_array_length(obj);
    if (arrLen != 3)
    {
        warn("Got %d items for vector '%s' in '%s', expected 3\n",
             arrLen,
             p->name,
             pname);
        return TRUE;
    }

    for (i = 0; i < 3; ++i)
    {
        tmp = (json_object*) array_list_get_idx(arr, i);
        if (!json_object_is_number(tmp))
        {
            warn("Got unexpected type '%s' in position %d "
                 "of key '%s' in '%s', expected number.\n",
                 showNBodyType(json_object_get_type_safe(tmp)),
                 i,
                 p->name,
                 pname);
            return TRUE;
        }

        ((real*) p->param)[i] = json_object_get_double(tmp);
    }

    return FALSE;
}

static bool readArray(const Parameter* p, const char* pname, json_object* obj, bool useDflt)
{
    int i, arrLen;
    array_list* arr;
    json_object* tmp;
    char* readArr;
    char* readLoc;
    NbodyReadFunc reader;

    reader = (NbodyReadFunc) p->conv;
    if (!reader || p->size == 0)
    {
        warn("Read function or size not set for array '%s' in '%s'\n", p->name, pname);
        return TRUE;
    }

    if (!json_object_is_type(obj, json_type_array))
    {
        typeWarning(nbody_type_array, p->name, pname, obj);
        return TRUE;
    }

    arr = json_object_get_array(obj);
    arrLen = json_object_array_length(obj);

    readArr = (char*) callocSafe(arrLen, p->size);

    for (i = 0; i < arrLen; ++i)
    {
        tmp = (json_object*) array_list_get_idx(arr, i);

        readLoc = readArr + i * p->size;  /* Index into mystery sized type  */
        if (reader((void*) readLoc, p->name, tmp))
        {
            warn("Failed to read array item in position %d "
                 "of key '%s' in '%s'.\n", i, p->name, pname);
            free(readArr);
            return TRUE;
        }
    }

    *((char**) p->param) = (char*) readArr;

    if (p->length)
        *p->length = arrLen;

    return FALSE;
}

static bool readEnum(const Parameter* p, const char* pname, json_object* obj, bool useDflt)
{
    generic_enum_t conv;
    ReadEnum reader;

    /* This is actually a json_type_string, which we read
     * into an enum, or take a default value */
    if (useDflt)
    {
        *((generic_enum_t*) p->param) = *((generic_enum_t*) p->dflt);
        return FALSE;
    }

    reader = (ReadEnum) p->conv;
    if (!reader)
    {
        warn("Error: read function not set for enum '%s'\n", p->name);
        return TRUE;
    }

    conv = reader(json_object_get_string(obj));
    if (conv == -1)
        return TRUE;

    *((generic_enum_t*) p->param) = conv;

    return FALSE;
}

static bool readObject(const Parameter* p, const char* pname, json_object* obj, bool useDflt)
{
    NbodyReadFunc readFunc = (NbodyReadFunc) p->conv;

    return readFunc(p->param, p->name, obj);
}

static bool readItem(const Parameter* p, const char* pname, json_object* obj, bool useDflt)
{
    switch (p->type)
    {
        case nbody_type_double:
            return readDouble(p, pname, obj, useDflt);

        case nbody_type_int:
            return readInt(p, pname, obj, useDflt);

        case nbody_type_boolean:
            return readBool(p, pname, obj, useDflt);

        case nbody_type_string:
            return readString(p, pname, obj, useDflt);

        case nbody_type_vector:
            return readVector(p, pname, obj, useDflt);

        case nbody_type_array:
            return readArray(p, pname, obj, useDflt);

        case nbody_type_object:
            return readObject(p, pname, obj, useDflt);

        case nbody_type_enum:
            return readEnum(p, pname, obj, useDflt);

        default:
            warn("Unhandled parameter type %d for key '%s' in '%s'\n",
                 p->type, p->name, pname);
            return 1;
    }
}

static void printAvailableFields(const Parameter* p, const char* pname)
{
    warn("Failed to read required item of correct type in group '%s'\n"
         "\tFields are:\n", pname);

    while (p->name)
    {
        warn("\t\t%s (%s)", p->name, showNBodyType(p->type));

        if (p->dflt)
            warn("  (optional)");
        warn(",\n");
        ++p;
    }
}

/* Read a set of related parameters, e.g. the main NBodyCtx.
   If the unique flag is set, it's only valid to have one of the items in the group.
   Otherwise, it tries to use all of the parameters, and warns when extra elements are found.
   Fails if it fails to find a parameter that isn't defaultable (returns nonzero).
 */
int readParameterGroup(const Parameter* g,   /* The set of parameters */
                       json_object* hdr,     /* The object of the group */
                       const char* pname)    /* parent name */
{
    const Parameter* p;
    const Parameter* q;
    json_object* obj;
    bool useDflt, defaultable = FALSE;
    bool found = FALSE, readError = FALSE;

    assert(g);
    assert(hdr);

    if (checkIsObject(hdr, pname))
        return TRUE;

    p = g;
    while (p->name && !readError)
    {
        defaultable = (p->dflt != NULL);
        useDflt = found = FALSE;

        obj = json_object_object_get(hdr, p->name);
        if (!obj)
        {
            if (defaultable)
                useDflt = TRUE;
            else
            {
                warn("Failed to find or got 'null' for required key '%s' in '%s'\n",
                     p->name,
                     pname);
                readError = TRUE;
                break;  /* abandon the loop and print useful debugging */
            }
        }
        else
        {
            found = TRUE;
        }

        readError = readItem(p, pname, obj, useDflt);

        /* Explicitly delete it so we can check for extra stuff */
        json_object_object_del(hdr, p->name);
        ++p;
    }

    /* Skip the extra parameter warning if there was an error, since
     * abandoning the loop leaves lots of stuff in it */
    if (!readError)
        warnExtraParams(hdr, pname);

    /* Report what was expected in more detail */
    if (readError || (!found && !defaultable))
    {
        printAvailableFields(g, pname);
        return 1;
    }

    return 0;
}

bool readTypedGroup(const ParameterSet* ps,    /* Set of possible options */
                    json_object* obj,          /* Object to read from */
                    const char* name,          /* Name of the item */
                    generic_enum_t* typeRead)  /* Read enum value */
{
    json_object* typeObj;
    nbody_type readType;
    const ParameterSet* p;

    if (checkIsObject(obj, name))
        return TRUE;

    typeObj = json_object_object_get(obj, "type");
    if (!typeObj)
    {
        warn("Didn't find type field in '%s'\n", name);
        return TRUE;
    }

    p = readParameterSet(ps, json_object_get_string(typeObj), name);

    /* Delete type field to prevent superfluous unknown field warning */
    json_object_object_del(obj, "type");
    if (!p)
        return TRUE;

    *typeRead = p->enumVal;
    return readParameterGroup(p->parameters, obj, name);
}

