/* ************************************************************************** */
/* CODE.H: define various extern things for code.c and io.c. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _CODE_H_
#define _CODE_H_

#include "defs.h"

typedef struct
{
    char* infile;        /* file name for snapshot input */
    char* outfile;       /* file name for snapshot output */
    real freq;           /* inverse of integration timestep */
    real freqout;        /* output frequency */
    real tstop;          /* time to stop calculation */
    char* headline;      /* message describing calculation */
    real tnow;           /* current value of time */
    real tout;           /* time of next output */
    int nstep;           /* number of time-steps */
    int nfcalc;          /* count force calculations */
    int n2bcalc;         /* count body-body interactions */
    int nbccalc;         /* count body-cell interactions */
    int nbody;           /* number of bodies in system */
    bodyptr bodytab;     /* points to array of bodies */
} NBodyCtx;

extern NBodyCtx ctx;

/* Global function prototypes. */

void initoutput(void);          /* open files for output */
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

