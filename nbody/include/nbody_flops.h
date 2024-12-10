#ifndef _NBODY_FLOPS_H_
#define _NBODY_FLOPS_H_

#include "nbody_types.h"

#ifdef NBODY_PAPI
#include <papi.h>

typedef struct {
    long long flops;
    long long mflops;
    double elapsed_time;
    int EventSet;
} NBFlopCounter;

/* Initialize PAPI counters */
int nbInitFlops(NBFlopCounter* counter);

/* Start counting flops */
int nbStartFlops(NBFlopCounter* counter);

/* Stop counting and get results */
int nbStopFlops(NBFlopCounter* counter);

/* Print FLOPS results */
void nbPrintFlops(const NBFlopCounter* counter);

#else

typedef struct {
    int dummy;  /* Empty struct when PAPI disabled */
} NBFlopCounter;

/* Stub functions when PAPI is disabled */
static inline int nbInitFlops(NBFlopCounter* counter) { return 0; }
static inline int nbStartFlops(NBFlopCounter* counter) { return 0; }
static inline int nbStopFlops(NBFlopCounter* counter) { return 0; }
static inline void nbPrintFlops(const NBFlopCounter* counter) {}

#endif /* NBODY_PAPI */

#endif /* _NBODY_FLOPS_H_ */
