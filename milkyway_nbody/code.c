/* ************************************************************************** */
/* CODE.C: hierarchical N-body code. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "code.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include <popt.h>
#include <json/json.h>

NBodyParams ps = { 0, };
Tree t = { 0, } ;
NBodyCtx ctx = { 0, };

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

    nbPms     = json_object_object_get_safe(hdr, "nbody-parameters");
    iniCoords = json_object_object_get_safe(hdr, "initial-coordinates");
    nbodyCtx  = json_object_object_get_safe(hdr, "nbody-context");
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
    ps.eps        = json_object_take_double(nbPms, "eps");
    ps.theta      = json_object_take_double(nbPms, "theta");


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


    /* Parameters related to the context */

    ctx.outfilename = json_object_take_string(nbodyCtx, "outfile");

    /* The json object has ownership of the string, so we need to copy it */
    initoutput(&ctx);

    ctx.headline    = json_object_take_string(nbodyCtx, "headline");
    ctx.tstop       = json_object_take_double(nbodyCtx, "tstop");
    ctx.nbody       = json_object_take_double(nbodyCtx, "nbody");
    ctx.allowIncest = json_object_take_bool(nbodyCtx, "allow-incest");
    ctx.usequad     = json_object_take_bool(nbodyCtx, "usequad");
    ctx.seed        = json_object_take_int(nbodyCtx, "seed");
    ctx.freq        = json_object_take_double(nbodyCtx, "freq");
    ctx.dtout       = json_object_take_double(nbodyCtx, "dtout");
    ctx.freqout     = json_object_take_double(nbodyCtx, "freqout");

    modelStr = json_object_take_string(nbodyCtx, "model");

    if (!strcasecmp(modelStr, "bh86"))
        ctx.model = BH86;
    else if (!strcasecmp(modelStr, "sw93"))
        ctx.model = SW93;
    else
    {
        fprintf(stderr, "Invalid model %s: Model options are either 'bh86' or 'sw93'\n", modelStr);
        exit(EXIT_FAILURE);
    }
    free(modelStr);

    /* use dt = dtnbody/2 to make sure we get enough orbit precision,
       this puts the results in ps.XC, ps.YC, ps.ZC, ps.XC, ps.YC, ps.ZC,
       which is then used by the testdata routine to do the shift.
       Be aware: changing the mass changes the orbit results, but this is OK */

    ps.eps = ps.r0 / (10 * sqrt((real)ctx.nbody));
    dtnbody = (1 / 10.0) * (1 / 10.0) * sqrt(((4 / 3) * M_PI * ps.r0 * ps.r0 * ps.r0) / (ps.PluMass));
    //dtnbody = pow(2.718,log(0.5)*kmax);
    ctx.freq = 1.0 / dtnbody;
    ctx.freqout = ctx.freq;
    ps.dtorbit = dtnbody / 2.0;


    /* Tree related parameters */
    t.rsize = json_object_take_double(treePms, "rsize");


    /* Scan through for leftover / unknown keys and provide warnings if any exist */
    warn_extra_params(nbodyCtx, "nbody-context");
    json_object_object_del(hdr, "nbody-context");

    warn_extra_params(nbPms, "nbody-parameters");
    json_object_object_del(hdr, "nbody-parameters");

    warn_extra_params(iniCoords, "initial-coordinates");
    json_object_object_del(hdr, "initial-coordinates");

    warn_extra_params(treePms, "tree");
    json_object_object_del(hdr, "tree");

    /* Now warn for entire groups on the header and whole file */
    warn_extra_params(hdr, "nbody-parameters-file");

    /* deref the top level object should take care of freeing whatever's left */
    json_object_put(obj);
}

static void initNBody(int argc, const char** argv)
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


char* headline = "Hierarchical N-body Code";   /* default id for run */

static void startrun(NBodyCtx*);           /* initialize system state */
static void stepsystem(void);         /* advance by one time-step */
static void testdata(void);           /* generate test data */
static void pickshell(vector, real);  /* pick point on shell */

/* MAIN: toplevel routine for hierarchical N-body code. */

int main(int argc, char* argv[])
{
    float chisqans = 0.0;

    initNBody(argc, (const char**) argv);

    // Calculate the reverse orbit
    printf("Calculating reverse orbit...");
    integrate();
    printf("done\n");

    printf("Beginning run...\n");
    startrun(&ctx);                 /* set params, input data */

    while (ctx.tnow < ctx.tstop - 1.0 / (1024 * ctx.freq)) /* while not past ctx.tstop */
        stepsystem();               /* advance N-body system */

    // Get the likelihood
    chisqans = chisq();
    printf("Run finished. chisq = %f\n", chisqans);

    stopoutput();               /* finish up output */

    free(ctx.headline);

    return 0;
}

/* STARTRUN: startup hierarchical N-body code. */

static void startrun(NBodyCtx* ctx)
{
    if (ctx->nbody < 1)              /* check input value */
        error("startrun: ctx.nbody = %d is absurd\n", ctx->nbody);

    //srand48((long) time(NULL));   /* set random generator */
    srand48((long) 0.0);    /* set random generator */
    testdata();             /* make test model */

    ctx->nstep = 0;                   /* start counting stps.eps */
    ctx->tout = ctx->tnow;            /* schedule first output */
}

/* TESTDATA: generate Plummer model initial conditions for test runs,
 * scaled to units such that M = -4E = G = 1 (Henon, Hegge, etc).
 * See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
 */

#define MFRAC  0.999                /* cut off 1-MFRAC of mass */

static void testdata(void)
{
    printf("Initializing plummer model...");
    real rsc, vsc, r, v, x, y;
    vector cmr, cmv;
    vector rshift, vshift, scaledrshift, scaledvshift;

    // The coordinates to shift the plummer sphere by
    (rshift)[0] = ps.XC;
    (rshift)[1] = ps.YC;
    (rshift)[2] = ps.ZC;
    (vshift)[0] = ps.XC;
    (vshift)[1] = ps.YC;
    (vshift)[2] = ps.ZC;

    printf("Shifting plummer sphere to r = (%f, %f, %f) v = (%f, %f, %f)...\n",
           rshift[0],
           rshift[1],
           rshift[2],
           vshift[0],
           vshift[1],
           vshift[2]);

    bodyptr p;

    ctx.headline = "Hierarchical code: Plummer model";
    /* supply default ctx.headline */
    ctx.tnow = 0.0;                 /* reset elapsed model time */
    ctx.bodytab = (bodyptr) allocate(ctx.nbody * sizeof(body));
    /* alloc space for bodies */
    rsc = ps.r0;               /* set length scale factor */
    vsc = rsqrt(ps.PluMass / rsc);         /* and recip. speed scale */
    CLRV(cmr);                  /* init cm pos, vel */
    CLRV(cmv);
    CLRV(scaledrshift);
    CLRV(scaledvshift);
    MULVS(scaledrshift, rshift, rsc);   /* Multiply shift by scale factor */
    MULVS(scaledvshift, vshift, vsc);   /* Multiply shift by scale factor */

    for (p = ctx.bodytab; p < ctx.bodytab + ctx.nbody; p++) /* loop over particles */
    {
        Type(p) = BODY;             /* tag as a body */
        Mass(p) = ps.PluMass / (real)ctx.nbody;            /* set masses equal */
        r = 1 / rsqrt(rpow(xrandom(0.0, MFRAC), /* pick r in struct units */
                           -2.0 / 3.0) - 1);
        pickshell(Pos(p), rsc * r);     /* pick scaled position */
        ADDV(Pos(p), Pos(p), rshift);       /* move the position */
        ADDV(cmr, cmr, Pos(p));         /* add to running sum */
        do                      /* select from fn g(x) */
        {
            x = xrandom(0.0, 1.0);      /* for x in range 0:1 */
            y = xrandom(0.0, 0.1);      /* max of g(x) is 0.092 */
        }
        while (y > x * x * rpow(1 - x * x, 3.5)); /* using von Neumann tech */
        v = rsqrt(2.0) * x / rpow(1 + r * r, 0.25); /* find v in struct units */
        pickshell(Vel(p), vsc * v);     /* pick scaled velocity */
        ADDV(Vel(p), Vel(p), vshift);       /* move the velocity */
        ADDV(cmv, cmv, Vel(p));         /* add to running sum */
    }
    DIVVS(cmr, cmr, (real) ctx.nbody);      /* normalize cm coords */
    DIVVS(cmv, cmv, (real) ctx.nbody);

    printf("done\n");
}

/* PICKSHELL: pick a random point on a sphere of specified radius. */

static void pickshell(vector vec, real rad)
{
    int k;
    real rsq, rsc;

    do                      /* pick point in NDIM-space */
    {
        for (k = 0; k < NDIM; k++)      /* loop over dimensions */
            vec[k] = xrandom(-1.0, 1.0);        /* pick from unit cube */
        DOTVP(rsq, vec, vec);           /* compute radius squared */
    }
    while (rsq > 1.0);                      /* reject if outside sphere */
    rsc = rad / rsqrt(rsq);         /* compute scaling factor */
    MULVS(vec, vec, rsc);           /* rescale to radius given */
}

/* STEPSYSTEM: advance N-body system one time-step. */

static void stepsystem(void)
{
    bodyptr p;
    real dt;
    vector dvel, dpos;
    int nfcalc;          /* count force calculations */
    int n2bcalc;         /* count body-body interactions */
    int nbccalc;         /* count body-cell interactions */


    if (ctx.nstep == 0)                 /* about to take 1st step? */
    {
        printf("Building tree...Starting Nbody simulation...\n");
        maketree(ctx.bodytab, ctx.nbody);       /* build tree structure */
        nfcalc = n2bcalc = nbccalc = 0;     /* zero counters */
        for (p = ctx.bodytab; p < ctx.bodytab + ctx.nbody; p++)
        {
            /* loop over all bodies */
            hackgrav(p, Mass(p) > 0.0);     /* get force on each */
            nfcalc++;               /* count force calcs */
            n2bcalc += ps.n2bterm;         /* and 2-body terms */
            nbccalc += ps.nbcterm;         /* and body-cell terms */
        }
        output();               /* do initial output */
    }
    dt = 1.0 / ctx.freq;                /* set basic time-step */
    for (p = ctx.bodytab; p < ctx.bodytab + ctx.nbody; p++) /* loop over all bodies */
    {
        MULVS(dvel, Acc(p), 0.5 * dt);      /* get velocity increment */
        ADDV(Vel(p), Vel(p), dvel);     /* advance v by 1/2 step */
        MULVS(dpos, Vel(p), dt);        /* get positon increment */
        ADDV(Pos(p), Pos(p), dpos);     /* advance r by 1 step */
    }

    maketree(ctx.bodytab, ctx.nbody);           /* build tree structure */
    nfcalc = n2bcalc = nbccalc = 0;     /* zero counters */
    for (p = ctx.bodytab; p < ctx.bodytab + ctx.nbody; p++) /* loop over bodies */
    {
        hackgrav(p, Mass(p) > 0.0);     /* get force on each */
        nfcalc++;               /* count force calcs */
        n2bcalc += ps.n2bterm;         /* and 2-body terms */
        nbccalc += ps.nbcterm;         /* and body-cell terms */
    }
    for (p = ctx.bodytab; p < ctx.bodytab + ctx.nbody; p++) /* loop over all bodies */
    {
        MULVS(dvel, Acc(p), 0.5 * dt);          /* get velocity increment */
        ADDV(Vel(p), Vel(p), dvel);             /* advance v by 1/2 step */
    }
    ctx.nstep++;                    /* count another time step */
    ctx.tnow = ctx.tnow + dt;           /* finally, advance time */
    output();                   /* do major or minor output */
}

