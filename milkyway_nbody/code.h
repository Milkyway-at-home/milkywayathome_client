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

void initoutput(NBodyCtx*);     /* open files for output */
void destroyCtx(NBodyCtx*);     /* close output files */
void inputdata(void);           /* read initial data file */
void maketree(bodyptr, int);    /* construct tree structure */
void hackgrav(bodyptr, bool);   /* compute force on body */
void output(void);              /* perform output operation */

/* Utility routines used in code.c and io.c.  These are defined in util.c
 * and getparam.c, which must be compiled with same choice of precision.
 */

real xrandom(real, real);         /* generate a random number */
void initparam(char**, char**);   /* initialize parameter pkg */
char* getparam(char*);            /* get parameter as char* */
int getiparam(char*);             /* get parameter as integer */
bool getbparam(char*);            /* get parameter as bool */
real getrparam(char*);            /* get parameter as real */

#endif /* _CODE_H_ */

