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

#include "json_params.h"
#include "code.h"

static void warn_extra_params(json_object* obj, const char* objName)
{
    json_object_object_foreach(obj,key,val)
    {
        fprintf(stderr,
                "Warning: In group '%s': Unknown field '%s': '%s'\n",
                objName,
                key,
                json_object_to_json_string(val));
    }
}

static inline json_object* json_object_object_get_safe(json_object* obj, const char* key)
{
    json_object* tmp;

    tmp = json_object_object_get(obj, key);
    if (!tmp)
    {
        fprintf(stderr, "Failed to find expected key '%s'\n", key);
        exit(EXIT_FAILURE);
    }
        return tmp;
}

/* Read named double with error checking from a json_object.
   We don't want to allow the automatic conversion to 0 for invalid objects.
   We also want to be accept integers as doubles.

   Also explicitly delete the object. Reference counting is unhappy with what I want to do.
   TODO: Maybe use the hash table directly to avoid a second lookup on the delete.
 */
static inline double json_object_take_double(json_object* obj, const char* key)
{
    double val;
    json_object* tmp;

    tmp = json_object_object_get_safe(obj, key);

    if (   json_object_is_type(tmp, json_type_double)
        || json_object_is_type(tmp, json_type_int))  /* It's OK to forget the decimal */
    {
        val = json_object_get_double(tmp);
        json_object_object_del(obj, key);
        return val;
    }

    fprintf(stderr, "Failed to read expected double for key '%s'\n", key);
    exit(EXIT_FAILURE);
}

/* Same as for json_object_take double.
   The string is copied, and needs to be freed.
 */
inline static char* json_object_take_string(json_object* obj, const char* key)
{
    json_object* tmp;
    char* str;

    tmp = json_object_object_get_safe(obj, key);
    if (json_object_is_type(tmp, json_type_null))
        str = NULL;
    else if (json_object_is_type(tmp, json_type_string))
    {
        /* The json_object owns the string, so we need to copy it */
        str = json_object_get_string(tmp);
        if (str)
            str = strdup(str);
    }
    else
    {
        fprintf(stderr, "Expected string or null for key '%s'\n", key);
        exit(EXIT_FAILURE);
    }

    json_object_object_del(obj, key);

    return str;
}

/* There is a pattern here. */
static inline bool json_object_take_bool(json_object* obj, const char* key)
{
    bool val;
    json_object* tmp;

    tmp = json_object_object_get_safe(obj, key);
    val = json_object_get_boolean(tmp);
    json_object_object_del(obj, key);

    return val;
}

static inline int json_object_take_int(json_object* obj, const char* key)
{
    int val;
    json_object* tmp;

    tmp = json_object_object_get_safe(obj, key);
    if (!json_object_is_type(tmp, json_type_int))
    {
        fprintf(stderr, "Expected int for key '%s'\n", key);
        exit(EXIT_FAILURE);
    }

    val = json_object_get_int(tmp);
    json_object_object_del(obj, key);
    return val;
}

static void get_params_from_json(json_object* obj)
{
    bool useRad = FALSE;
    char* modelStr;
    static real dtnbody;
    static real kmax;

    json_object* tmp;
    json_object* nbPms;
    json_object* hdr;
    json_object* iniCoords;
    json_object* nbodyCtx;
    json_object* treePms;

    /* Check that this is actually one of our files */
    if (   !json_object_is_type(obj, json_type_object)
        || !(hdr = json_object_object_get(obj, "nbody-parameters-file")))
    {
        fprintf(stderr, "Parameters not in expected format.\n");
        exit(EXIT_FAILURE);
    }

    nbodyCtx  = json_object_object_get_safe(hdr, "nbody-context");
    nbPms     = json_object_object_get_safe(hdr, "nbody-parameters");
    iniCoords = json_object_object_get_safe(hdr, "initial-coordinates");
    treePms   = json_object_object_get_safe(hdr, "tree");

    if (!(nbPms && iniCoords && nbodyCtx && treePms))
    {
        fprintf(stderr, "Parameter group missing\n");
        exit(EXIT_FAILURE);
    }


    /* First group of parameters */
    ps.PluMass    = json_object_take_double(nbPms, "PluMass");
    ps.r0         = json_object_take_double(nbPms, "r0");
    ps.orbittstop = json_object_take_double(nbPms, "orbittstop");


    /* Parameters related to initial coordinates */
    useRad = json_object_take_bool(iniCoords, "angle-use-radians");

    ps.lstart = json_object_take_double(iniCoords, "lstart");
    ps.bstart = json_object_take_double(iniCoords, "bstart");

    if (useRad)
    {
        ps.lstart = d2r(ps.lstart);
        ps.bstart = d2r(ps.bstart);
    }

    ps.Rstart    = json_object_take_double(iniCoords, "Rstart");
    ps.Xinit     = json_object_take_double(iniCoords, "Xinit");
    ps.Yinit     = json_object_take_double(iniCoords, "Yinit");
    ps.Zinit     = json_object_take_double(iniCoords, "Zinit");
    ps.sunGCDist = json_object_take_double(iniCoords, "sunGCDist");

    ps.Xinit = ps.Rstart * cos(ps.lstart) * cos(ps.bstart) - ps.sunGCDist;
    ps.Yinit = ps.Rstart * sin(ps.lstart) * cos(ps.bstart);
    ps.Zinit = ps.Rstart * sin(ps.bstart);


    /* TODO: Maybe seed */
    /* Parameters related to the context */

    ctx.headline    = json_object_take_string(nbodyCtx, "headline");
    ctx.outfilename = json_object_take_string(nbodyCtx, "outfile");
    ctx.nbody       = json_object_take_double(nbodyCtx, "nbody");
    ctx.seed        = json_object_take_int(nbodyCtx, "seed");
    ctx.tstop       = json_object_take_double(nbodyCtx, "tstop");
    ctx.allowIncest = json_object_take_bool(nbodyCtx, "allow-incest");
    ctx.usequad     = json_object_take_bool(nbodyCtx, "use-quadruople-corrections");
    ctx.freq        = json_object_take_double(nbodyCtx, "freq");
    ctx.dtout       = json_object_take_double(nbodyCtx, "dtout");
    ctx.freqout     = json_object_take_double(nbodyCtx, "freqout");
    ctx.theta       = json_object_take_double(nbodyCtx, "accuracy-parameter");
    ctx.eps         = json_object_take_double(nbodyCtx, "potential-softening");

    /* The json object has ownership of the string, so we need to copy it */
    initoutput(&ctx);
    modelStr = json_object_take_string(nbodyCtx, "criterion");
    if (!strcasecmp(modelStr, "bh86"))
        ctx.criterion = BH86;
    else if (!strcasecmp(modelStr, "sw93"))
        ctx.criterion = SW93;
    else
    {
        fprintf(stderr, "Invalid model %s: Model options are either 'bh86' or 'sw93'\n", modelStr);
        exit(EXIT_FAILURE);
    }
    free(modelStr);

    /* Scan through for leftover / unknown keys and provide warnings if any exist */
    warn_extra_params(nbodyCtx, "nbody-context");
    json_object_object_del(hdr, "nbody-context");



    /* use dt = dtnbody/2 to make sure we get enough orbit precision,
       this puts the results in ps.XC, ps.YC, ps.ZC, ps.XC, ps.YC, ps.ZC,
       which is then used by the testdata routine to do the shift.
       Be aware: changing the mass changes the orbit results, but this is OK */

    ctx.eps = ps.r0 / (10 * sqrt((real)ctx.nbody));
    dtnbody = (1 / 10.0) * (1 / 10.0) * sqrt(((4 / 3) * M_PI * ps.r0 * ps.r0 * ps.r0) / (ps.PluMass));
    //dtnbody = pow(2.718,log(0.5)*kmax);
    ctx.freq = 1.0 / dtnbody;
    ctx.freqout = ctx.freq;
    ps.dtorbit = dtnbody / 2.0;


    /* Tree related parameters */
    t.rsize = json_object_take_double(treePms, "rsize");
    warn_extra_params(treePms, "tree");
    json_object_object_del(hdr, "tree");


    warn_extra_params(nbPms, "nbody-parameters");
    json_object_object_del(hdr, "nbody-parameters");

    warn_extra_params(iniCoords, "initial-coordinates");
    json_object_object_del(hdr, "initial-coordinates");

    /* Now warn for entire groups on the header and whole file */
    warn_extra_params(hdr, "nbody-parameters-file");

    /* deref the top level object should take care of freeing whatever's left */
    json_object_put(obj);
}

void initNBody(int argc, const char** argv)
{
    poptContext context;
    int o;
    static char* inputFile = NULL;  /* input JSON file */
    static char* inputStr = NULL;   /* a string of JSON to use directly */
    static json_object* obj;

    static const struct poptOption options[] =
    {
        {
            "input-file", 'f',
            POPT_ARG_STRING, &inputFile,
            0, "Input file to read", NULL
        },

        {
            "input-string", 's',
            POPT_ARG_STRING, &inputStr,
            0, "Input given as string", NULL
        },

        POPT_AUTOHELP

        { NULL, 0, 0, NULL, 0, NULL, NULL }
    };


    context = poptGetContext(argv[0],
                             argc,
                             argv,
                             options,
                             POPT_CONTEXT_POSIXMEHARDER);

    if (argc < 2)
    {
        poptPrintUsage(context, stderr, 0);
        poptFreeContext(context);
        exit(EXIT_FAILURE);
    }

    while ( ( o = poptGetNextOpt(context)) >= 0 );

    /* Check for invalid options, and must have one of input file or input string */
    if (     o < -1
         ||  (inputFile && inputStr)
         || !(inputFile || inputStr)
        )
    {
        poptPrintHelp(context, stderr, 0);
        exit(EXIT_FAILURE);
    }

    poptFreeContext(context);


    if (inputFile)
    {
        obj = json_object_from_file(inputFile);
        if (is_error(obj))
        {
            fprintf(stderr,
                    "Failed to read file '%s'. Perhaps not found, or a parse error?\n",
                    inputFile);
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        obj = json_tokener_parse(inputStr);
        if (is_error(obj))
        {
            fprintf(stderr, "Failed to parse given string\n");
            exit(EXIT_FAILURE);
        }

    }

    get_params_from_json(obj);
}


