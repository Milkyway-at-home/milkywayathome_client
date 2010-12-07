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

#define histogramPhi 128.79
#define histogramTheta 54.39
#define histogramPsi 90.70
#define histogramStartRaw ((real) -50.0)
#define histogramEndRaw ((real) 50.0)
#define histogramBinSize ((real) 2.9411764705882355)
#define histogramCenter ((real) 0.0)


/* static const real nanN = NAN; */
static real nanN;

static int readParameterGroup(const Parameter* g, json_object* hdr, const char* parentName);


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

static nbody_type json_object_get_type_safe(json_object* obj)
{
    return obj ? json_object_get_type(obj) : nbody_type_null;
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

static disk_t readDiskT(const char* str)
{
    if (!strcasecmp(str, "miyamoto-nagai"))
        return MiyamotoNagaiDisk;
    if (!strcasecmp(str, "exponential"))
        return ExponentialDisk;
    else
        warn("Invalid disk type %s: Disk options are either 'miyamoto-nagai', "
             "or 'exponential'\n",
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

static halo_t readHaloT(const char* str)
{
    if (!strcasecmp(str, "logarithmic"))
        return DwarfModelPlummer;
    if (!strcasecmp(str, "nfw"))
        return DwarfModelKing;
    if (!strcasecmp(str, "triaxial"))
        return DwarfModelDehnen;
    else
        warn("Invalid model %s: Model options are "
             "'logarithmic', 'nfw', or 'triaxial'\n", str);

    return -1;
}

/* From a given set of options with different parameters, find the matching set */
static const ParameterSet* readParameterSet(const ParameterSet* ps, const char* str, const char* name)
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


static bool readTypedGroup(const ParameterSet* ps,    /* Set of possible options */
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

static bool readInitialConditions(InitialConditions* ic, const char* pname, json_object* obj)
{
    const InitialConditions defaultIC =
        {
            /* .position     */ EMPTY_VECTOR,
            /* .velocity     */ EMPTY_VECTOR,
            /* .useGalC      */ FALSE,
            /* .useRadians   */ FALSE,
            /* .reverseOrbit */ TRUE
        };

    const Parameter initialConditionParams[] =
        {
            BOOL_PARAM("use-galactic-coordinates", &ic->useGalC),
            BOOL_PARAM_DFLT("angle-use-radians", &ic->useRadians, &defaultIC.useRadians),
            BOOL_PARAM_DFLT("reverse-orbit", &ic->reverseOrbit, &defaultIC.reverseOrbit),
            VEC_PARAM("position", &ic->position),
            VEC_PARAM("velocity", &ic->velocity),
            NULLPARAMETER
        };

    return readParameterGroup(initialConditionParams, obj, pname);
}

static bool readHaloParams(Halo* halo, const char* pname, json_object* obj)
{
    const Parameter nfwParams[] =
        {
            DBL_PARAM("scale-length", &halo->scale_length),
            DBL_PARAM("vhalo",        &halo->vhalo),
            NULLPARAMETER
        };

    const Parameter logarithmicParams[] =
        {
            DBL_PARAM("vhalo",        &halo->vhalo),
            DBL_PARAM("z-flattening", &halo->flattenZ),
            DBL_PARAM("scale-length", &halo->scale_length),
            NULLPARAMETER
        };

    const Parameter triaxialParams[] =
        {
            DBL_PARAM("vhalo",          &halo->vhalo),
            DBL_PARAM("scale-length",   &halo->scale_length),
            DBL_PARAM("x-flattening",   &halo->flattenX),
            DBL_PARAM("y-flattening",   &halo->flattenY),
            DBL_PARAM("z-flattening",   &halo->flattenZ),
            DBL_PARAM("triaxial-angle", &halo->triaxAngle),
            NULLPARAMETER
        };

    const ParameterSet haloOptions[] =
        {
            { "logarithmic", LogarithmicHalo, logarithmicParams },
            { "nfw",         NFWHalo,         nfwParams },
            { "triaxial",    TriaxialHalo,    triaxialParams },
            NULL_PARAMETER_SET
        };

    return readTypedGroup(haloOptions, obj, pname, (generic_enum_t*) &halo->type);
}

static bool readDwarfModel(DwarfModel* model, const char* parentName, json_object* obj)
{
    const bool defaultIgnore = FALSE;

    /* The current different dwarf models all use the same parameters */
    const Parameter dwarfModelParams[] =
        {
            /* FIXME: Hack: Defaulting on NAN's so we can ignore them
             * in the file, to be filled in by the server sent
             * FitParams. This will probably result in unfortunate
             * things when using the file. */
            ENUM_PARAM("type",               &model->type,           (NbodyReadFunc) readDwarfModelT),
            INT_PARAM("nbody",               &model->nbody),
            DBL_PARAM_DFLT("mass",           &model->mass,           &nanN),
            DBL_PARAM_DFLT("scale-radius",   &model->scale_radius,   &nanN),
            DBL_PARAM_DFLT("timestep",       &model->timestep,       &nanN),
            DBL_PARAM_DFLT("orbit-timestep", &model->orbit_timestep, &nanN),
            BOOL_PARAM_DFLT("ignore-final",  &model->ignoreFinal, &defaultIgnore),
            OBJ_PARAM("initial-conditions",  &model->initialConditions, (NbodyReadFunc) readInitialConditions),
            NULLPARAMETER
        };

    return readParameterGroup(dwarfModelParams, obj, parentName);
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

static bool json_object_is_number(json_object* obj)
{
    return obj && (json_object_is_type(obj, json_type_double) || json_object_is_type(obj, json_type_int));
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
        warn("Error: expected number for '%s' in '%s', but got %s\n",
             p->name,
             pname,
             showNBodyType(json_object_get_type_safe(obj)));
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
             showNBodyType(json_object_get_type_safe(obj)));
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
             showNBodyType(json_object_get_type_safe(obj)));
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
             showNBodyType(json_object_get_type_safe(obj)));
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
        warn("Error: expected type vector for '%s' in '%s', but got %s\n",
             p->name,
             pname,
             showNBodyType(json_object_get_type_safe(obj)));
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
        warn("Error: expected type array for '%s' in '%s', but got %s\n",
             p->name,
             pname,
             showNBodyType(json_object_get_type_safe(obj)));
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

static bool readDiskParams(Disk* disk, const char* pname, json_object* obj)
{
    const Parameter miyamotoParams[] =
        {
            DBL_PARAM("mass",         &disk->mass),
            DBL_PARAM("scale-length", &disk->scale_length),
            DBL_PARAM("scale-height", &disk->scale_height),
            NULLPARAMETER
        };

    const Parameter exponentialParams[] =
        {
            DBL_PARAM("mass",         &disk->mass),
            DBL_PARAM("scale-length", &disk->scale_length),
            NULLPARAMETER
        };

    const ParameterSet diskOptions[] =
        {
            { "exponential",    ExponentialDisk,   exponentialParams },
            { "miyamoto-nagai", MiyamotoNagaiDisk, miyamotoParams },
            NULL_PARAMETER_SET
        };

    return readTypedGroup(diskOptions, obj, pname, (generic_enum_t*) &disk->type);
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
static int readParameterGroup(const Parameter* g,         /* The set of parameters */
                              json_object* hdr,           /* The object of the group */
                              const char* parentName)     /* non-null if within another group */
{
    const Parameter* p;
    const Parameter* q;
    bool useDflt, defaultable = FALSE;
    json_object* obj;
    const char* pname;
    bool found = FALSE, done = FALSE, readError = FALSE;
    int subRc;
    nbody_type readType;
    NbodyReadFunc readFunc;

    /* Name of the group we're in */
    pname = parentName ? parentName : "<root>";

    assert(g);
    assert(hdr);

    p = g;

    if (checkIsObject(hdr, pname))
        return TRUE;

    /* CHECKME: Handling of defaultable, cases in unique? */
    while (p->name && !done && !readError)
    {
        defaultable = (p->dflt != NULL);
        useDflt = FALSE;
        found = FALSE;

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

            case nbody_type_object:
                readFunc = (NbodyReadFunc) p->conv;
                readError = readFunc(p->param, p->name, obj);
                break;

            case nbody_type_enum:
                readError = readEnum(p, pname, obj, useDflt);
                break;

            default:
                warn("Unhandled parameter type %d for key '%s' in '%s'\n",
                     p->type, p->name, pname);
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

    /* Report what was expected in more detail */
    if (   (!found && !defaultable)
        || readError)
    {
        warn("Failed to find required item of correct type in group '%s'\n"
             "\tFields are:\n", pname);
        q = g;
        while (q->name)
        {
            warn("\t\t%s (%s)", q->name, showNBodyType(q->type));

            if (q->dflt)
                warn("  (optional)");
            warn(",\n");
            ++q;
        }

        return 1;
    }

    return 0;
}

static bool readSphericalParams(Spherical* sph, const char* pname, json_object* obj)
{
    const Parameter sphericalParams[] =
        {
            DBL_PARAM("mass",     &sph->mass),
            DBL_PARAM("r0-scale", &sph->scale),
            NULLPARAMETER
        };

    const ParameterSet sphericalOptions[] =
        {
            { "sphere", SphericalPotential, sphericalParams },
            NULL_PARAMETER_SET
        };

    return readTypedGroup(sphericalOptions, obj, pname, (generic_enum_t*) &sph->type);
}

static bool readPotential(Potential* pot, const char* pname, json_object* obj)
{
    const Parameter potentialItems[] =
        {
            OBJ_PARAM("disk",      &pot->disk,      (NbodyReadFunc) readDiskParams),
            OBJ_PARAM("halo",      &pot->halo,      (NbodyReadFunc) readHaloParams),
            OBJ_PARAM("spherical", &pot->sphere[0], (NbodyReadFunc) readSphericalParams),
            NULLPARAMETER
        };

    return readParameterGroup(potentialItems, obj, pname);
}

static bool readNbodyContext(NBodyCtx* ctx, const char* pname, json_object* obj)
{
    /* Constants used for defaulting. Each field only used if
     * specified in the actual parameter tables. */
    const NBodyCtx defaultCtx =
        {
            /* Grr lack of C99 named struct initializers in MSVC */
            /* .pot             */  EMPTY_POTENTIAL,
            /* .models          */  NULL,
            /* .modelNum        */  0,
            /* .nbody           */  0,
            /* .timestep        */  0.0,
            /* .time_evolve     */  0.0,
            /* .orbit_timestep  */  0.0,
            /* time_orbit       */  0.0,

            /* .headline        */  NULL,
            /* .outfilename     */  NULL,
            /* .histogram       */  NULL,
            /* .histout         */  NULL,
            /* .outfile         */  NULL,

            /* .freqout         */  0.0,
            /* .theta           */  0.0,
            /* .eps2            */  0.0,

            /* .tree_rsize      */  4.0,
            /* .sunGCDist       */  8.0,
            /* .criterion       */  NEWCRITERION,
            /* .seed            */  0,
            /* .usequad         */  TRUE,
            /* .allowIncest     */  FALSE,

            /* .outputCartesian */  FALSE,
            /* .outputBodies    */  FALSE,
            /* .outputHistogram */  FALSE,
            /* .cp_filename     */  NULL,
            /* .cp_resolved     */  ""
        };

    /* Must be null terminated arrays */
    const Parameter nbodyCtxParams[] =
        {
            STR_PARAM("headline",                    &ctx->headline),
            INT_PARAM_DFLT("seed",                   &ctx->seed, &defaultCtx.seed),
            BOOL_PARAM("use-quadrupole-corrections", &ctx->usequad),

            BOOL_PARAM_DFLT("allow-incest",          &ctx->allowIncest, &defaultCtx.allowIncest),
            DBL_PARAM_DFLT("accuracy-parameter",     &ctx->theta, &defaultCtx.theta),
            ENUM_PARAM_DFLT("criterion",             &ctx->criterion, &defaultCtx.criterion, (NbodyReadFunc) readCriterion),
            DBL_PARAM_DFLT("eps2", &ctx->eps2, &nanN),

            OBJ_PARAM("potential", &ctx->pot, (NbodyReadFunc) readPotential),
            DBL_PARAM_DFLT("time-evolve",    &ctx->time_evolve,    &nanN),
            DBL_PARAM_DFLT("time-orbit",     &ctx->time_orbit,     &nanN),

            DBL_PARAM_DFLT("timestep",       &ctx->timestep,       &nanN),
            DBL_PARAM_DFLT("orbit_timestep", &ctx->orbit_timestep, &nanN),

            ARRAY_PARAM("dwarf-model", &ctx->models, sizeof(DwarfModel), &ctx->modelNum, (NbodyReadFunc) readDwarfModel),
            DBL_PARAM_DFLT("sun-gc-dist", &ctx->sunGCDist, &defaultCtx.sunGCDist),
            DBL_PARAM_DFLT("tree_rsize", &ctx->tree_rsize, &defaultCtx.tree_rsize),
            NULLPARAMETER
        };

    return readParameterGroup(nbodyCtxParams, obj, pname);
}

static bool readHistogramParams(HistogramParams* hist, const char* pname, json_object* obj)
{
    const HistogramParams defaultHistogram =
        {
            /* .phi      */  histogramPhi,
            /* .theta    */  histogramTheta,
            /* .psi      */  histogramPsi,
            /* .startRaw */  histogramStartRaw,
            /* .endRaw   */  histogramEndRaw,
            /* .binSize  */  histogramBinSize,
            /* .center   */  histogramCenter
        };

    const Parameter histogramParams[] =
        {
            DBL_PARAM_DFLT("phi",     &hist->phi,      &defaultHistogram.phi),
            DBL_PARAM_DFLT("theta",   &hist->theta,    &defaultHistogram.theta),
            DBL_PARAM_DFLT("psi",     &hist->psi,      &defaultHistogram.psi),
            DBL_PARAM_DFLT("start",   &hist->startRaw, &defaultHistogram.startRaw),
            DBL_PARAM_DFLT("end",     &hist->endRaw,   &defaultHistogram.endRaw),
            DBL_PARAM_DFLT("binsize", &hist->binSize,  &defaultHistogram.binSize),
            DBL_PARAM_DFLT("center",  &hist->center,   &defaultHistogram.center),
            NULLPARAMETER
        };

    *hist = defaultHistogram;    /* Set all items to default */

    return readParameterGroup(histogramParams, obj, pname);
}


/* CHECKME: Assumption that all enums are the same size, size of an
 * int, which I think is always true on all compilers, but I'm
 * probably wrong. */

/* TODO: Could use hash table more directly and avoid all the extra
 * lookups, but it probably doesn't matter. */


/* Read the parameters from the top level json object into ctx. It
 * destroys the object in the process. */
int getParamsFromJSON(NBodyCtx* ctx,         /* Context to fill */
                      HistogramParams* hist, /* Histogram parameters to fill */
                      json_object* fileObj)  /* Parsed JSON file */

{
    const Parameter parameters[] =
        {
            OBJ_PARAM("nbody-context", ctx,  (NbodyReadFunc) readNbodyContext),
            OBJ_PARAM("histogram",     hist, (NbodyReadFunc) readHistogramParams),
            NULLPARAMETER
        };

    json_object* hdr;
    int rc;

    nanN = NAN; /* Work around MSVC stupidity. Actually set value of nan that's defaulted to */

    /* Check that this is actually one of our files */
    if (   !json_object_is_type(fileObj, nbody_type_object)
        || !(hdr = json_object_object_get(fileObj, "nbody-parameters-file")))
    {
        warn("Parameters not in expected format.\n");
        return 1;
    }

    /* loop through table of accepted sets of parameters */
    rc = readParameterGroup(parameters, hdr, NULL);

    /* deref the top level object should take care of freeing whatever's left */
    json_object_put(fileObj);

    return rc;
}

