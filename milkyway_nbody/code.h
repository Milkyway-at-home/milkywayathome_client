/* ************************************************************************** */
/* CODE.H: define various extern things for code.c and io.c. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _CODE_H_
#define _CODE_H_

#include "defs.h"

#define OUTFILENAMELEN 1024

typedef struct
{
    char* headline;      /* message describing calculation */
    model_t model;       /* bh86 or sw93 */
    bool usequad;        /* use quadrupole corrections */
    bool allowIncest;
    real freq;           /* inverse of integration timestep */
    real freqout;        /* output frequency */
    real tstop;          /* time to stop calculation */
    real tnow;           /* current value of time */
    real tout;           /* time of next output */
    int nbody;           /* number of bodies in system */
    int nstep;           /* number of time-steps */
    bodyptr bodytab;     /* points to array of bodies */
    FILE* outfile;                      /* file for snapshot output */
    char outfilename[OUTFILENAMELEN];   /* filename for snapshot output */
} NBodyCtx;

extern NBodyParams ps;
extern Tree t;
extern NBodyCtx ctx;

/* Global function prototypes. */

void initoutput(NBodyCtx*);     /* open files for output */
void stopoutput(void);          /* close output files */
void inputdata(void);           /* read initial data file */
void maketree(bodyptr, int);    /* construct tree structure */
void hackgrav(bodyptr, bool);   /* compute force on body */
void output(void);              /* perform output operation */

/*  * Utility routines used in code.c and io.c.  These are defined in util.c
 * and getparam.c, which must be compiled with same choice of precision.
 */

real xrandom(real, real);         /* generate a random number */
void initparam(char**, char**);   /* initialize parameter pkg */
char* getparam(char*);            /* get parameter as char* */
int getiparam(char*);             /* get parameter as integer */
bool getbparam(char*);            /* get parameter as bool */
real getrparam(char*);            /* get parameter as real */

#endif /* _CODE_H_ */

