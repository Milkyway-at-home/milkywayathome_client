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
#include "code.h"
#include "json_params.h"
#include "defs.h"


/* Read command line arguments and initialize the context and state */
void initNBody(int argc, const char** argv)
{
    poptContext context;
    int o;
    static char* inputFile = NULL;  /* input JSON file */
    static char* inputStr = NULL;   /* a string of JSON to use directly */
    static json_object* obj;

    /* FIXME: There's a small leak of the inputFile from use of
       poptGetNextOpt().  Some mailing list post suggestst that this
       is some kind of semi-intended bug to work around something or
       other while maintaining ABI compatability */
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
         || !(inputFile || inputStr))
    {
        poptPrintHelp(context, stderr, 0);
        free(inputFile);
        free(inputStr);
        poptFreeContext(context);
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
            free(inputFile);
            exit(EXIT_FAILURE);
        }
        free(inputFile);
    }
    else
    {
        obj = json_tokener_parse(inputStr);
        free(inputStr);
        if (is_error(obj))
            fail("Failed to parse given string\n");
    }

    get_params_from_json(obj);

}

/* Reads a name of the criterion into the C value */
static criterion_t readCriterion(const char* str)
{
    if (!strcasecmp(str, "bh86"))
        return BH86;
    else if (!strcasecmp(str, "sw93"))
        return SW93;
    else
        fail("Invalid model %s: Model options are either 'bh86' or 'sw93'\n", str);
}

/* Iterate through remaining keys in obj to provide useful warnings of
 * unknown parameters in the file */
static void warn_extra_params(json_object* obj, const char* grpName)
{
    json_object_object_foreach(obj,key,val)
    {
        fprintf(stderr,
                "Warning: In group '%s': Unknown field '%s': '%s'\n",
                grpName,
                key,
                json_object_to_json_string(val));
    }
}

/* Read a set of related parameters, e.g. the main NBodyCtx */
inline static void readParameterGroup(const ParameterGroup* g, json_object* hdr)
{
    const Parameter* p;
    bool defaultable, useDflt;
    json_object* grp;
    json_object* obj;

    p = g->parameters;

    grp = json_object_object_get(hdr, g->name);
    if (!grp)
        fail("Failed to find expected group '%s'\n", g->name);

    /* set of associated parameters */
    while (p->name)
    {
        useDflt = FALSE;
        defaultable = (p->dflt != NULL);

        obj = json_object_object_get(grp, p->name);
        if (!obj)   /* not found */
        {
            if (defaultable)
                useDflt = TRUE;
            else
            {
                fail("Failed to find required key '%s' in group '%s'\n",
                     p->name,
                     g->name);
            }
        }
        else   /* found the object */
        {
            if (json_object_is_type(obj, json_type_null))
            {
                if (defaultable)
                    useDflt = TRUE;
                else
                {
                    /* Found it, but trying to give null when a value is required */
                    fail("Got null for required key '%s' in group '%s'\n",
                         p->name,
                         g->name);
                }
            }
        }

        /* TODO: Type checking in need of much improvement */
        switch (p->type)
        {
            case json_type_int:
                if (useDflt)
                {
                    *((int*) p->param) = *((int*) p->dflt);
                    break;
                }

                /* If we expect an int, and we find a string, it's an enum */
                if (p->conv)
                {
                    if (!json_object_is_type(obj, json_type_string))
                        fail("Expected enum, did not get string "
                             "for key '%s' in group '%s'\n",
                             p->name,
                             g->name);

                    *((int*) p->param) = (int) p->conv(json_object_get_string(obj));
                    break;
                }

                *((int*) p->param) = json_object_get_int(obj);
                break;

            case json_type_double:
                /* Cast to real is important */
                *((real*) p->param) = useDflt ? *((double*) p->dflt) : json_object_get_double(obj);
                break;

            case json_type_string:
                /* json_object has ownership of string, so copy it */
                *((char**) p->param) =
                    useDflt ? *((char**) p->dflt) : strdup(json_object_get_string(obj));
                break;

            case json_type_boolean:
                *((bool*) p->param) = useDflt ? *((bool*) p->dflt) : json_object_get_boolean(obj);
                break;

            default:
                fail("Unhandled parameter type %d for key '%s' in group '%s'\n",
                     p->type,
                     p->name,
                     g->name);
        }

        /* Explicitly delete this from the json_object so we can check
         * for unknown/extra parameters. Deref is not enough here. */
        json_object_object_del(grp, p->name);
        ++p;
    }

    warn_extra_params(grp, g->name);

}

/* CHECKME: Assumption that all enums are the same size, size of an
 * int, which I think is always true on all compilers, but I'm
 * probably wrong. */

/* TODO: Could use hash table more directly and avoid all the extra
 * lookups, but it probably doesn't matter. */

/* Read the parameters from the top level json object. It destroys the
 * object in the process. */
void get_params_from_json(json_object* fileObj)
{
    /* Constants used for defaulting. Each field only used if
     * specified in the actual parameter tables. */
    static const NBodyCtx defaultCtx =
    {
        .pot = { { 0 }, {0,0,0,0} , {0,0,0,0,0,0,0}, NULL },
        .nbody = 1024,
        .tstop = NAN,
        .dtout = NAN,
        .freq = NAN,
        .freqout = NAN,
        .usequad = TRUE,
        .allowIncest = FALSE,
        .criterion = BH86,
        .theta = 0.0,
        .eps = 0.0,
        .seed = 0,
        .outfile = NULL,
        .outfilename = NULL,
        .headline = NULL
    };

    static const Tree defaultTree =
    {
        .rsize = 4.0
    };

    /* Must be null terminated arrays */
    static const Parameter nbodyCtxParams[] =
        {
            { "outfile", json_type_string, &ctx.outfilename, NULL, NULL },
            { "headline", json_type_string, &ctx.headline, NULL, NULL },
            { "accuracy-parameter", json_type_double, &ctx.theta, &defaultCtx.theta, NULL },
            { "criterion", json_type_int, &ctx.criterion, NULL, readCriterion },
            { "use-quadrupole-corrections", json_type_boolean, &ctx.usequad, NULL, NULL },
            { "allow-incest", json_type_boolean, &ctx.allowIncest, NULL, NULL },

            { "seed", json_type_int, &ctx.seed, NULL, NULL },
            { "nbody", json_type_int, &ctx.nbody, NULL, NULL },
            { NULL, json_type_null, NULL, NULL, NULL }
        };

    static const Parameter treeParams[] =
    {
        { "rsize", json_type_double, &t.rsize, &defaultTree.rsize, NULL }
    };

    static const ParameterGroup parameters[] =
    {
        { "nbody-context", nbodyCtxParams },
        { "tree", treeParams },
        // { "initial-conditions", { 0 } },
        { NULL, { NULL, json_type_null, NULL, NULL, NULL } }
    };

    ParameterGroup* g;
    json_object* hdr;

    /* Check that this is actually one of our files */
    if (   !json_object_is_type(fileObj, json_type_object)
        || !(hdr = json_object_object_get(fileObj, "nbody-parameters-file")))
    {
        fprintf(stderr, "Parameters not in expected format.\n");
        exit(EXIT_FAILURE);
    }

    /* loop through table of accepted sets of parameters */
    g = parameters;
    while (g->name)
    {
        readParameterGroup(g, hdr);

        /* Explicitly delete the group so we can check for extra stuff */
        json_object_object_del(hdr, g->name);
        ++g;
    }

    /* Warnings for unknown groups */
    warn_extra_params(hdr, "nbody-parameters-file");

    printf("End read parameters\n");
    printContext(&ctx);

    /* deref the top level object should take care of freeing whatever's left */
    json_object_put(fileObj);

    //json_object_object_del(fileObj, "nbody-parameters-file");

}


