/* ************************************************************************** */
/* code.h: define various extern things for main.c and io.c. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _NBODY_PRIV_H_
#define _NBODY_PRIV_H_

#define _GNU_SOURCE

#include <json/json.h>

#include "nbody_config.h"
#include "nbody_types.h"
#include "json_params.h"
#include "stdinc.h"
#include "vectmath.h"
#include "real.h"
#include "nbody_util.h"

#if BOINC_APPLICATION
  #include <boinc_api.h>
  #if BOINC_DEBUG
    #include <diagnostics.h>
  #endif /* BOINC_DEBUG */
#endif /* BOINC_APPLICATION */

/* Global function prototypes. */

void initoutput(NBodyCtx*);             /* open files for output */
void openCheckpoint(NBodyCtx* ctx);
void thawState(const NBodyCtx* ctx, NBodyState* st);
void nbody_ctx_destroy(NBodyCtx* ctx);  /* close output files */
void nbody_state_destroy(NBodyState* st);

void inputdata(void);                   /* read initial data file */
void maketree(const NBodyCtx*, NBodyState*);    /* construct tree structure */
void hackgrav(const NBodyCtx*, NBodyState*, bodyptr, bool);   /* compute force on body */
void output(const NBodyCtx* ctx, const NBodyState* st);  /* perform output operation */
void boincOutput(const NBodyCtx* ctx, const NBodyState* st);
void nbody_boinc_output(const NBodyCtx* ctx, const NBodyState* st);
void generatePlummer(const NBodyCtx* ctx, const InitialConditions* ic, NBodyState* st);

void integrate(const NBodyCtx* ctx, InitialConditions* ic);

void acceleration(vectorptr restrict acc, const NBodyCtx* ctx, const vectorptr restrict pos);

void sphericalAccel(vectorptr restrict acc, const Spherical* sph, const vectorptr restrict pos);
void miyamotoNagaiDiskAccel(vectorptr restrict acc, const Disk* d, const vectorptr restrict pos);
void exponentialDiskAccel(vectorptr restrict acc, const Disk* d, const vectorptr restrict pos);
void triaxialHaloAccel(vectorptr restrict acc, const Halo* h, const vectorptr restrict pos);
void logHaloAccel(vectorptr restrict acc, const Halo* h, const vectorptr restrict pos);
void nfwHaloAccel(vectorptr restrict acc, const Halo* h, const vectorptr restrict pos);

/* Utility routines used in code.c and io.c.  These are defined in util.c
 * and getparam.c, which must be compiled with same choice of precision.
 */

real xrandom(real, real);         /* generate a random number */
void initparam(char**, char**);   /* initialize parameter pkg */
char* getparam(char*);            /* get parameter as char* */
int getiparam(char*);             /* get parameter as integer */
bool getbparam(char*);            /* get parameter as bool */
real getrparam(char*);            /* get parameter as real */

/* Types -> String */
const char* showBool(const bool);
const char* showCriterionT(const criterion_t);
const char* showSphericalT(const spherical_t);
const char* showDiskT(const disk_t);
const char* showHaloT(const halo_t);
const char* showDwarfModelT(const dwarf_model_t);
char* showSpherical(const Spherical*);
char* showHalo(const Halo*);
char* showDisk(const Disk*);
char* showPotential(const Potential*);
char* showDwarfModel(const DwarfModel*);
char* showInitialConditions(const InitialConditions*);
char* showContext(const NBodyCtx*);
char* showVector(const vector v);
void printContext(const NBodyCtx*);
void printInitialConditions(const InitialConditions*);

#endif /* _NBODY_PRIV_H_ */

