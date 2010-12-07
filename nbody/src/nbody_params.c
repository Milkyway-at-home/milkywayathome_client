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
#include "nbody_params.h"
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

    const MWParameter initialConditionParams[] =
        {
            BOOL_PARAM("use-galactic-coordinates", &ic->useGalC),
            BOOL_PARAM_DFLT("angle-use-radians", &ic->useRadians, &defaultIC.useRadians),
            BOOL_PARAM_DFLT("reverse-orbit", &ic->reverseOrbit, &defaultIC.reverseOrbit),
            VEC_PARAM("position", &ic->position),
            VEC_PARAM("velocity", &ic->velocity),
            NULL_MWPARAMETER
        };

    return mwReadParameterGroup(initialConditionParams, obj, pname);
}

static bool readHaloParams(Halo* halo, const char* pname, json_object* obj)
{
    const MWParameter nfwParams[] =
        {
            DBL_PARAM("scale-length", &halo->scale_length),
            DBL_PARAM("vhalo",        &halo->vhalo),
            NULL_MWPARAMETER
        };

    const MWParameter logarithmicParams[] =
        {
            DBL_PARAM("vhalo",        &halo->vhalo),
            DBL_PARAM("z-flattening", &halo->flattenZ),
            DBL_PARAM("scale-length", &halo->scale_length),
            NULL_MWPARAMETER
        };

    const MWParameter triaxialParams[] =
        {
            DBL_PARAM("vhalo",          &halo->vhalo),
            DBL_PARAM("scale-length",   &halo->scale_length),
            DBL_PARAM("x-flattening",   &halo->flattenX),
            DBL_PARAM("y-flattening",   &halo->flattenY),
            DBL_PARAM("z-flattening",   &halo->flattenZ),
            DBL_PARAM("triaxial-angle", &halo->triaxAngle),
            NULL_MWPARAMETER
        };

    const MWParameterSet haloOptions[] =
        {
            { "logarithmic", LogarithmicHalo, logarithmicParams },
            { "nfw",         NFWHalo,         nfwParams },
            { "triaxial",    TriaxialHalo,    triaxialParams },
            NULL_MWPARAMETERSET
        };

    return mwReadTypedGroup(haloOptions, obj, pname, (generic_enum_t*) &halo->type);
}

static bool readDwarfModel(DwarfModel* model, const char* parentName, json_object* obj)
{
    const bool defaultIgnore = FALSE;

    /* The current different dwarf models all use the same parameters */
    const MWParameter dwarfModelParams[] =
        {
            /* FIXME: Hack: Defaulting on NAN's so we can ignore them
             * in the file, to be filled in by the server sent
             * FitParams. This will probably result in unfortunate
             * things when using the file. */
            ENUM_PARAM("type",               &model->type,           (MWReadFunc) readDwarfModelT),
            INT_PARAM("nbody",               &model->nbody),
            DBL_PARAM_DFLT("mass",           &model->mass,           &nanN),
            DBL_PARAM_DFLT("scale-radius",   &model->scale_radius,   &nanN),
            DBL_PARAM_DFLT("timestep",       &model->timestep,       &nanN),
            DBL_PARAM_DFLT("orbit-timestep", &model->orbit_timestep, &nanN),
            BOOL_PARAM_DFLT("ignore-final",  &model->ignoreFinal, &defaultIgnore),
            OBJ_PARAM("initial-conditions",  &model->initialConditions, (MWReadFunc) readInitialConditions),
            NULL_MWPARAMETER
        };

    return mwReadParameterGroup(dwarfModelParams, obj, parentName);
}

static bool readDiskParams(Disk* disk, const char* pname, json_object* obj)
{
    const MWParameter miyamotoParams[] =
        {
            DBL_PARAM("mass",         &disk->mass),
            DBL_PARAM("scale-length", &disk->scale_length),
            DBL_PARAM("scale-height", &disk->scale_height),
            NULL_MWPARAMETER
        };

    const MWParameter exponentialParams[] =
        {
            DBL_PARAM("mass",         &disk->mass),
            DBL_PARAM("scale-length", &disk->scale_length),
            NULL_MWPARAMETER
        };

    const MWParameterSet diskOptions[] =
        {
            { "exponential",    ExponentialDisk,   exponentialParams },
            { "miyamoto-nagai", MiyamotoNagaiDisk, miyamotoParams },
            NULL_MWPARAMETERSET
        };

    return mwReadTypedGroup(diskOptions, obj, pname, (generic_enum_t*) &disk->type);
}

static bool readSphericalParams(Spherical* sph, const char* pname, json_object* obj)
{
    const MWParameter sphericalParams[] =
        {
            DBL_PARAM("mass",     &sph->mass),
            DBL_PARAM("r0-scale", &sph->scale),
            NULL_MWPARAMETER
        };

    const MWParameterSet sphericalOptions[] =
        {
            { "sphere", SphericalPotential, sphericalParams },
            NULL_MWPARAMETERSET
        };

    return mwReadTypedGroup(sphericalOptions, obj, pname, (generic_enum_t*) &sph->type);
}

static bool readPotential(Potential* pot, const char* pname, json_object* obj)
{
    const MWParameter potentialItems[] =
        {
            OBJ_PARAM("disk",      &pot->disk,      (MWReadFunc) readDiskParams),
            OBJ_PARAM("halo",      &pot->halo,      (MWReadFunc) readHaloParams),
            OBJ_PARAM("spherical", &pot->sphere[0], (MWReadFunc) readSphericalParams),
            NULL_MWPARAMETER
        };

    return mwReadParameterGroup(potentialItems, obj, pname);
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
    const MWParameter nbodyCtxParams[] =
        {
            STR_PARAM("headline",                    &ctx->headline),
            INT_PARAM_DFLT("seed",                   &ctx->seed, &defaultCtx.seed),
            BOOL_PARAM("use-quadrupole-corrections", &ctx->usequad),

            BOOL_PARAM_DFLT("allow-incest",          &ctx->allowIncest, &defaultCtx.allowIncest),
            DBL_PARAM_DFLT("accuracy-parameter",     &ctx->theta, &defaultCtx.theta),
            ENUM_PARAM_DFLT("criterion",             &ctx->criterion, &defaultCtx.criterion, (MWReadFunc) readCriterion),
            DBL_PARAM_DFLT("eps2", &ctx->eps2, &nanN),

            OBJ_PARAM("potential", &ctx->pot, (MWReadFunc) readPotential),
            DBL_PARAM_DFLT("time-evolve",    &ctx->time_evolve,    &nanN),
            DBL_PARAM_DFLT("time-orbit",     &ctx->time_orbit,     &nanN),

            DBL_PARAM_DFLT("timestep",       &ctx->timestep,       &nanN),
            DBL_PARAM_DFLT("orbit_timestep", &ctx->orbit_timestep, &nanN),

            ARRAY_PARAM("dwarf-model", &ctx->models, sizeof(DwarfModel), &ctx->modelNum, (MWReadFunc) readDwarfModel),
            DBL_PARAM_DFLT("sun-gc-dist", &ctx->sunGCDist, &defaultCtx.sunGCDist),
            DBL_PARAM_DFLT("tree_rsize", &ctx->tree_rsize, &defaultCtx.tree_rsize),
            NULL_MWPARAMETER
        };

    return mwReadParameterGroup(nbodyCtxParams, obj, pname);
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

    const MWParameter histogramParams[] =
        {
            DBL_PARAM_DFLT("phi",     &hist->phi,      &defaultHistogram.phi),
            DBL_PARAM_DFLT("theta",   &hist->theta,    &defaultHistogram.theta),
            DBL_PARAM_DFLT("psi",     &hist->psi,      &defaultHistogram.psi),
            DBL_PARAM_DFLT("start",   &hist->startRaw, &defaultHistogram.startRaw),
            DBL_PARAM_DFLT("end",     &hist->endRaw,   &defaultHistogram.endRaw),
            DBL_PARAM_DFLT("binsize", &hist->binSize,  &defaultHistogram.binSize),
            DBL_PARAM_DFLT("center",  &hist->center,   &defaultHistogram.center),
            NULL_MWPARAMETER
        };

    *hist = defaultHistogram;    /* Set all items to default */

    return mwReadParameterGroup(histogramParams, obj, pname);
}


/* CHECKME: Assumption that all enums are the same size, size of an
 * int, which I think is always true on all compilers, but I'm
 * probably wrong. */

/* TODO: Could use hash table more directly and avoid all the extra
 * lookups, but it probably doesn't matter. */

/* FIXME: Duplicates in a group lead to bad things */

/* Read the parameters from the top level json object into ctx. It
 * destroys the object in the process. */
int nbodyGetParamsFromJSON(NBodyCtx* ctx,         /* Context to fill */
                           HistogramParams* hist, /* Histogram parameters to fill */
                           json_object* fileObj)  /* Parsed JSON file */

{
    const MWParameter parameters[] =
        {
            OBJ_PARAM("nbody-context", ctx,  (MWReadFunc) readNbodyContext),
            OBJ_PARAM("histogram",     hist, (MWReadFunc) readHistogramParams),
            NULL_MWPARAMETER
        };

    json_object* hdr;
    int rc;

    nanN = NAN; /* Work around MSVC stupidity. Actually set value of nan that's defaulted to */

    /* Check that this is actually one of our files */
    if (   !json_object_is_type(fileObj, mw_type_object)
        || !(hdr = json_object_object_get(fileObj, "nbody-parameters-file")))
    {
        warn("Parameters not in expected format.\n");
        return 1;
    }

    /* loop through table of accepted sets of parameters */
    rc = mwReadParameterGroup(parameters, hdr, "<root>");

    /* deref the top level object should take care of freeing whatever's left */
    json_object_put(fileObj);

    return rc;
}

