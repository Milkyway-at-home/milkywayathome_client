#include "nbody_flops.h"
#include "milkyway_util.h"

#ifdef NBODY_PAPI

int nbInitFlops(NBFlopCounter* counter)
{
    int retval;
    
    /* Initialize PAPI library */
    if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT)
    {
        mw_printf("PAPI library initialization error: %d\n", retval);
        return 1;
    }

    /* Create event set */
    counter->EventSet = PAPI_NULL;
    if ((retval = PAPI_create_eventset(&counter->EventSet)) != PAPI_OK)
    {
        mw_printf("PAPI create eventset error: %d\n", retval);
        return 1;
    }

    /* Add FLOPS counter */
    if ((retval = PAPI_add_event(counter->EventSet, PAPI_FP_OPS)) != PAPI_OK)
    {
        mw_printf("PAPI add FLOPS event error: %d\n", retval);
        return 1;
    }

    counter->flops = 0;
    counter->mflops = 0;
    counter->elapsed_time = 0.0;

    return 0;
}

int nbStartFlops(NBFlopCounter* counter)
{
    int retval;
    if ((retval = PAPI_start(counter->EventSet)) != PAPI_OK)
    {
        mw_printf("PAPI start error: %d\n", retval);
        return 1;
    }
    return 0;
}

int nbStopFlops(NBFlopCounter* counter)
{
    int retval;
    long long values[1];

    if ((retval = PAPI_stop(counter->EventSet, values)) != PAPI_OK)
    {
        mw_printf("PAPI stop error: %d\n", retval);
        return 1;
    }

    counter->flops = values[0];
    counter->elapsed_time = PAPI_get_real_usec() / 1e6;
    counter->mflops = (long long)(counter->flops / (counter->elapsed_time * 1e6));

    return 0;
}

void nbPrintFlops(const NBFlopCounter* counter)
{
    mw_printf("Total FLOPS: %lld\n", counter->flops);
    mw_printf("MFLOPS: %lld\n", counter->mflops);
    mw_printf("Time: %.3f seconds\n", counter->elapsed_time);
}

#endif /* NBODY_PAPI */