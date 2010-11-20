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

static const real nanN = NAN;

static int readParameterGroup(const Parameter* g, json_object* hdr, const Parameter* parent, generic_enum_t* group_type);


/* also works for json_object_type */
static const char* showNBodyType(nbody_type bt)
{
    static const char* table[] =
        {
            [nbody_type_null]       = "null",
            [nbody_type_boolean]    = "bool",
            [nbody_type_double]     = "double",
            [nbody_type_int]        = "int",
            [nbody_type_object]     = "object",
            [nbody_type_array]      = "array",
            [nbody_type_vector]     = "vector",
            [nbody_type_string]     = "string",
            [nbody_type_enum]       = "enum",
            [nbody_type_group]      = "group",
            [nbody_type_group_item] = "group_item"
        };

    if (bt > nbody_type_group_item)
    {
        warn("Trying to show unknown nbody_type %d\n", bt);
        return NULL;
    }

    return table[bt];
}

static bool readGroupItem(const Parameter* p, const char* pname, json_object* obj, generic_enum_t* group_type)
{
    if (!p->dflt)
    {
        warn("Expected nbody_type_group_item for "
             "'%s' in '%s', but no enum value set\n", p->name, pname);
        return TRUE;
    }

    if (!group_type)
    {
        warn("Group does not have a type (Trying to read the root?)\n");
        return TRUE;
    }

    *group_type = *((generic_enum_t*) p->dflt);

    return FALSE;
}

/* Reads a name of the criterion into the C value, with name str */
static criterion_t readCriterion(const char* str)
{
    if (!strcasecmp(str, "new-criterion"))
        return NEWCRITERION;
    if (!strcasecmp(str, "exact"))
        return EXACT;
    if (!strcasecmp(str, "bh86"))
        return BH86;
    else if (!strcasecmp(str, "sw93"))
        return SW93;
    else
        warn("Invalid criterion %s: Criterion options are either 'bh86', "
             "'sw93', 'exact' or 'new-criterion' (default),\n",
             str);

    return -1;
}

static dwarf_model_t readDwarfModelT(const char* str)
{
    if (!strcasecmp(str, "plummer"))
        return DwarfModelPlummer;
    if (!strcasecmp(str, "king"))
        return DwarfModelKing;
    if (!strcasecmp(str, "dehnen"))
        return DwarfModelDehnen;
    else
        warn("Invalid model %s: Model options are "
             "'plummer', 'king', 'dehnen'\n", str);

    return -1;
}

static bool readDwarfModel(DwarfModel* model, json_object* obj, const Parameter* parent)
{
    const InitialConditions defaultIC =
        {
            .useGalC      = FALSE,
            .useRadians   = FALSE,
            .reverseOrbit = TRUE,
            .position     = EMPTY_VECTOR,
            .velocity     = EMPTY_VECTOR
        };

    const Parameter initialConditionParams[] =
        {
            BOOL_PARAM("useGalC", &model->initialConditions.useGalC),
            BOOL_PARAM_DFLT("angle-use-radians", &model->initialConditions.useRadians, &defaultIC.useRadians),
            BOOL_PARAM_DFLT("reverse-orbit", &model->initialConditions.reverseOrbit, &defaultIC.reverseOrbit),
            VEC_PARAM("position", &model->initialConditions.position),
            VEC_PARAM("velocity", &model->initialConditions.velocity),
            NULLPARAMETER
        };

    /* The current different dwarf models all use the same parameters */
    const Parameter dwarfModelParams[] =
        {
            /* FIXME: Hack: Defaulting on NAN's so we can ignore them
             * in the file, to be filled in by the server sent
             * FitParams. This will probably result in unfortunate
             * things when using the file. */
            ENUM_PARAM("type",               &model->type,           (GenericReadFunc) readDwarfModelT),
            INT_PARAM("nbody",               &model->nbody),
            DBL_PARAM_DFLT("mass",           &model->mass,           &nanN),
            DBL_PARAM_DFLT("scale-radius",   &model->scale_radius,   &nanN),
            DBL_PARAM_DFLT("timestep",       &model->timestep,       &nanN),
            DBL_PARAM_DFLT("orbit-timestep", &model->orbit_timestep, &nanN),
            OBJ_PARAM("initial-conditions",  initialConditionParams),
            NULLPARAMETER
        };

    return readParameterGroup(dwarfModelParams, obj, parent, NULL);
}

static void printParameter(Parameter* p)
{
    if (p)
    {
        printf("parameter = {\n"
               "  name  = %s\n"
               "  type  = %d\n"
               "  param = %p\n"
               "  dflt  = %p\n"
               "  conv  = %p\n"
               "  uniq  = %d\n"
               "  pms   = %p\n"
               "};\n",
               p->name,
               p->type,
               p->param,
               p->dflt,
               p->conv,
               p->unique,
               p->parameters);
    }
    else
    {
        printf("<NULL PARAMETER>\n");
    }
}


/* Iterate through remaining keys in obj to provide useful warnings of
 * unknown parameters in the file. Returns true if any found. */
static bool warnExtraParams(json_object* obj, const char* grpName)
{
    bool haveExtra = FALSE;

    json_object_object_foreach(obj,key,val)
    {
        haveExtra = TRUE;
        fprintf(stderr,
                "Warning: In group '%s': Unknown field '%s': '%s'\n",
                grpName,
                key,
                json_object_to_json_string(val));
    }

    return haveExtra;
}


static bool readDouble(const Parameter* p, const char* pname, json_object* obj, bool useDflt)
{
    /* json_type_int and double are OK for numbers. i.e. you can leave off the decimal.
       We don't want the other conversions, which just give you 0.0 for anything else.
    */

    if (useDflt)
        *((real*) p->param) = *((real*) p->dflt);
    else if (   json_object_is_type(obj, json_type_double)
             || json_object_is_type(obj, json_type_int))
    {
        *((real*) p->param) = (real) json_object_get_double(obj);
    }
    else
    {
        warn("Error: expected number for '%s' in '%s', but got %s\n",
             p->name,
             pname,
             showNBodyType(p->type));
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
        warn("Error: expected type int for '%s' in '%s', but got %s\n",
             p->name,
             pname,
             showNBodyType(json_object_get_type(obj)));
        return TRUE;
    }

    return FALSE;
}

static bool readBool(const Parameter* p, const char* pname, json_object* obj, bool useDflt)
{
    if (useDflt)
        *((bool*) p->param) = *((int*) p->dflt);
    else if (json_object_is_type(obj, json_type_boolean))
        *((bool*) p->param) = (bool) json_object_get_boolean(obj);
    else
    {
        warn("Error: expected type boolean for '%s' in '%s', but got %s\n",
             p->name,
             pname,
             showNBodyType(json_object_get_type(obj)));
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
        warn("Error: expected type string for '%s' in '%s', but got %s\n",
             p->name,
             pname,
             showNBodyType(json_object_get_type(obj)));
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
    assert(json_object_is_type(obj, json_type_array));
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

        if (   !json_object_is_type(tmp, json_type_double)
            && !json_object_is_type(tmp, json_type_int))
        {
            warn("Got unexpected type '%s' in position %d "
                 "of key '%s' in '%s', expected number.\n",
                 showNBodyType(json_object_get_type(tmp)),
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
    ArrayRead reader;

    reader = (ArrayRead) p->conv;
    if (!reader || p->size == 0)
    {
        warn("Read function or size not set for array '%s' in '%s'\n", p->name, pname);
        return TRUE;
    }

    assert(json_object_is_type(obj, json_type_array));

    arr = json_object_get_array(obj);
    arrLen = json_object_array_length(obj);

    readArr = (char*) callocSafe(arrLen, p->size);

    for (i = 0; i < arrLen; ++i)
    {
        tmp = (json_object*) array_list_get_idx(arr, i);

        readLoc = readArr + i * p->size;  /* Index into mystery sized type  */
        if (reader((void*) readLoc, tmp, p))
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

/* Read a set of related parameters, e.g. the main NBodyCtx.
   If the unique flag is set, it's only valid to have one of the items in the group.
   Otherwise, it tries to use all of the parameters, and warns when extra elements are found.
   Fails if it fails to find a parameter that isn't defaultable (returns nonzero).
 */
static int readParameterGroup(const Parameter* g,        /* The set of parameters */
                              json_object* hdr,           /* The object of the group */
                              const Parameter* parent,    /* non-null if within another group */
                              generic_enum_t* group_type) /* non-null if within a group. i.e. type of model */
{
    const Parameter* p;
    const Parameter* q;
    bool useDflt, defaultable = FALSE;
    json_object* obj;
    const char* pname;
    bool unique;
    bool found = FALSE, done = FALSE, readError = FALSE;
    int subRc;

    if (parent)
    {
        pname = parent->name;     /* Name of the group we're in */
        unique = parent->unique;  /* Expect one, or try to take all of them*/
    }
    else
    {
        pname = "<root>";       /* Probably at the root of the configuration */
        unique = FALSE;
    }

    assert(g);
    assert(hdr);

    p = g;

    /* CHECKME: Handling of defaultable, cases in unique? */
    while (p->name && !done && !readError)
    {
        defaultable = (p->dflt != NULL);
        useDflt = FALSE;
        found = FALSE;

        obj = json_object_object_get(hdr, p->name);
        if (!obj)
        {
            if (unique)  /* We get another chance */
            {
                ++p;
                continue;
            }

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
            if (unique)    /* Only want one parameter, everything else is extra and a warning */
                done = TRUE;
        }

        /* TODO: Better type checking might be nice. Mostly now relies
         * on not screwing up the tables. */
        switch (p->type)
        {
            case nbody_type_double:
                readError = readDouble(p, pname, obj, useDflt);
                break;

            case nbody_type_int:
                readError = readInt(p, pname, obj, useDflt);
                break;

            case nbody_type_boolean:  /* CHECKME: Size */
                readError = readBool(p, pname, obj, useDflt);
                break;

            case nbody_type_string:
                readError = readString(p, pname, obj, useDflt);
                break;

            case nbody_type_vector:
                readError = readVector(p, pname, obj, useDflt);
                break;

            case nbody_type_array:
                readError = readArray(p, pname, obj, useDflt);
                break;

            case nbody_type_group_item:
                readError = readGroupItem(p, pname, obj, group_type);

                /* fall through to nbody_type_object */
            case nbody_type_group:
            case nbody_type_object:
                /* Fail now if sub group fails */
                if ((subRc = readParameterGroup(p->parameters, obj, p, (generic_enum_t*) p->param)))
                    return subRc;
                break;

            case nbody_type_enum:
                readError = readEnum(p, pname, obj, useDflt);
                break;

            default:
                warn("Unhandled parameter type %d for key '%s' in '%s'\n",
                     p->type,
                     p->name,
                     pname);
                return 1;
        }

        /* Explicitly delete it so we can check for extra stuff */
        json_object_object_del(hdr, p->name);

        ++p;
    }

    /* Skip the extra parameter warning if there was an error, since
     * abandoning the loop leaves lots of stuff in it */
    if (!readError)
        warnExtraParams(hdr, pname);

    /* FIXME: This condition is confusing and probably could be better */
    /* Report what was expected in more detail */
    if (   ((!found || readError) && !defaultable)
        || (!found && unique) )
    {
        warn("Failed to find required item of correct type in group '%s'\n", pname);

        if (unique)
            warn("\tExpected to find one of the following:\n");
        else
            warn("\tFields are:\n");

        q = g;
        while (q->name)
        {
            warn("\t\t%s (%s)", q->name, showNBodyType(q->type));

            if (q->dflt && !unique)
                warn("  (optional)");
            warn(",\n");
            ++q;
        }

        return 1;
    }

    return 0;
}

/* CHECKME: Assumption that all enums are the same size, size of an
 * int, which I think is always true on all compilers, but I'm
 * probably wrong. */

/* TODO: Could use hash table more directly and avoid all the extra
 * lookups, but it probably doesn't matter. */


/* Read the parameters from the top level json object into ctx. It
 * destroys the object in the process. */
int getParamsFromJSON(NBodyCtx* ctx,               /* Context to fill */
                      json_object* fileObj)        /* Parsed JSON file */

{
    /* Constants used for defaulting. Each field only used if
     * specified in the actual parameter tables. */
    const NBodyCtx defaultCtx =
        {
            .pot = EMPTY_POTENTIAL,
            .freqout = NAN,
            .usequad = TRUE,
            .sunGCDist  = 8.0,
            .allowIncest = FALSE,
            .outputCartesian = FALSE,
            .criterion = NEWCRITERION,
            .theta = 0.0,
            .seed = 0,
            .tree_rsize = 4.0,
            .outfile = NULL,
            .outfilename = NULL,
            .headline = NULL
        };

    /* Spherical potential options */
    const Parameter sphericalParams[] =
        {
            DBL_PARAM("mass",     &ctx->pot.sphere[0].mass),
            DBL_PARAM("r0-scale", &ctx->pot.sphere[0].scale),
            NULLPARAMETER
        };

    /* TODO: Different spherical potentials, more spherical potentials */
    /* This is maybe redundant and annoying, (spherical in spherical),
     * but to be consistent with the others. Also part of the way to
     * having multiple spherical potentials. */
    const spherical_t sphT = SphericalPotential;
    const Parameter sphericalOptions[] =
        {
            GROUP_PARAM_ITEM("sphere", &sphT,  sphericalParams),
            NULLPARAMETER
        };

    /* Disk potential options */
    const Parameter miyamotoParams[] =
        {
            DBL_PARAM("mass",         &ctx->pot.disk.mass),
            DBL_PARAM("scale-length", &ctx->pot.disk.scale_length),
            DBL_PARAM("scale-height", &ctx->pot.disk.scale_height),
            NULLPARAMETER
        };

    const Parameter exponentialParams[] =
        {
            DBL_PARAM("mass",         &ctx->pot.disk.mass),
            DBL_PARAM("scale-length", &ctx->pot.disk.scale_length),
            NULLPARAMETER
        };

    /* Can't take the address of a literal. Also using pointers to
     * guarantee a consistent size. */
    const disk_t mnd = MiyamotoNagaiDisk, expd = ExponentialDisk;
    const Parameter diskOptions[] =
        {
            GROUP_PARAM_ITEM("miyamoto-nagai", &mnd,  miyamotoParams),
            GROUP_PARAM_ITEM("exponential",    &expd, exponentialParams),
            NULLPARAMETER
        };

    /* Halo options */
    const Parameter nfwParams[] =
        {
            DBL_PARAM("scale-length", &ctx->pot.halo.scale_length),
            DBL_PARAM("vhalo",        &ctx->pot.halo.vhalo),
            NULLPARAMETER
        };

    const Parameter logarithmicParams[] =
        {
            DBL_PARAM("vhalo",        &ctx->pot.halo.vhalo),
            DBL_PARAM("z-flattening", &ctx->pot.halo.flattenZ),
            DBL_PARAM("scale-length", &ctx->pot.halo.scale_length),
            NULLPARAMETER
        };

    const Parameter triaxialParams[] =
        {
            DBL_PARAM("vhalo",          &ctx->pot.halo.vhalo),
            DBL_PARAM("scale-length",   &ctx->pot.halo.scale_length),
            DBL_PARAM("x-flattening",   &ctx->pot.halo.flattenX),
            DBL_PARAM("y-flattening",   &ctx->pot.halo.flattenY),
            DBL_PARAM("z-flattening",   &ctx->pot.halo.flattenZ),
            DBL_PARAM("triaxial-angle", &ctx->pot.halo.triaxAngle),
            NULLPARAMETER
        };

    const halo_t logHaloT = LogarithmicHalo, nfwHaloT = NFWHalo, triHaloT = TriaxialHalo;
    const Parameter haloOptions[] =
        {
            GROUP_PARAM_ITEM("nfw",         &nfwHaloT, nfwParams),
            GROUP_PARAM_ITEM("logarithmic", &logHaloT, logarithmicParams),
            GROUP_PARAM_ITEM("triaxial",    &triHaloT, triaxialParams),
            NULLPARAMETER
        };

    const Parameter potentialItems[] =
        {
            GROUP_PARAM("disk",      &ctx->pot.disk.type,      diskOptions),
            GROUP_PARAM("halo",      &ctx->pot.halo.type,      haloOptions),
            GROUP_PARAM("spherical", &ctx->pot.sphere[0].type, sphericalOptions),
            NULLPARAMETER
        };

    /* Must be null terminated arrays */
    const Parameter nbodyCtxParams[] =
        {
            STR_PARAM("headline",                    &ctx->headline),
            INT_PARAM_DFLT("seed",                   &ctx->seed, &defaultCtx.seed),
            BOOL_PARAM("use-quadrupole-corrections", &ctx->usequad),

            BOOL_PARAM_DFLT("allow-incest",          &ctx->allowIncest, &defaultCtx.allowIncest),
            DBL_PARAM_DFLT("accuracy-parameter",     &ctx->theta, &defaultCtx.theta),
            ENUM_PARAM_DFLT("criterion",             &ctx->criterion, &defaultCtx.criterion, (GenericReadFunc) readCriterion),
            DBL_PARAM_DFLT("eps2", &ctx->eps2, &nanN),

            OBJ_PARAM("potential", potentialItems),
            DBL_PARAM_DFLT("time-evolve",    &ctx->time_evolve,    &nanN),
            DBL_PARAM_DFLT("time-orbit",     &ctx->time_orbit,     &nanN),
            ARRAY_PARAM("dwarf-model", &ctx->models, sizeof(DwarfModel), &ctx->modelNum, (GenericReadFunc) readDwarfModel),
            DBL_PARAM_DFLT("sun-gc-dist", &ctx->sunGCDist, &defaultCtx.sunGCDist),
            DBL_PARAM_DFLT("tree_rsize", &ctx->tree_rsize, &defaultCtx.tree_rsize),
            NULLPARAMETER
        };

    const Parameter parameters[] =
        {
            OBJ_PARAM("nbody-context", nbodyCtxParams),
            NULLPARAMETER
        };

    json_object* hdr;
    int rc;

    /* Check that this is actually one of our files */
    if (   !json_object_is_type(fileObj, nbody_type_object)
        || !(hdr = json_object_object_get(fileObj, "nbody-parameters-file")))
    {
        warn("Parameters not in expected format.\n");
        return 1;
    }

    /* loop through table of accepted sets of parameters */
    rc = readParameterGroup(parameters, hdr, NULL, NULL);

    /* deref the top level object should take care of freeing whatever's left */
    json_object_put(fileObj);


    return rc;
}

