/* ************************************************************************** */
/* CODE.H: define various extern things for code.c and io.c. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _CODE_H_
#define _CODE_H_

#include "defs.h"

extern Tree t;
extern NBodyCtx ctx;
extern NBodyState st;
extern NBodyParams ps;

/* Global function prototypes. */

void initoutput(NBodyCtx*);             /* open files for output */
void nbody_ctx_destroy(NBodyCtx* ctx);  /* close output files */
void inputdata(void);                   /* read initial data file */
void maketree(bodyptr, int);    /* construct tree structure */
void hackgrav(bodyptr, bool);   /* compute force on body */
void output(void);              /* perform output operation */

void acceleration(real* posvec, real* accvec);
void integrate();

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
const char* showBool(bool);
const char* showCriterionT(criterion_t);
const char* showSphericalT(spherical_t);
const char* showDiskT(disk_t);
const char* showHaloT(halo_t);
const char* showDwarfModelT(dwarf_model_t);
char* showSpherical(Spherical*);
char* showHalo(Halo*);
char* showDisk(Disk*);
char* showPotential(Potential*);
char* showDwarfModel(DwarfModel*);
char* showInitialConditions(InitialConditions*);
char* showContext(NBodyCtx*);
void printContext(NBodyCtx*);
void printInitialConditions(InitialConditions*);

#endif /* _CODE_H_ */

