#ifndef TBX_TIMING_H
#define TBX_TIMING_H

#include <sys/types.h>
#include <sys/time.h>
#if defined(CRLIBM_TYPECPU_ITANIUM) && defined(__INTEL_COMPILER)
#include<ia64intrin.h>
#endif

/*
 * For the user interested in High resolution timer, to my knowledge 
 * the best library is "RTAI", unfortunately it works only on Linux OS
 */


typedef union u_tbx_tick
{
  unsigned long long tick;

  struct
  {
    unsigned long low;
    unsigned long high;
  } sub;

  struct timeval timev;
} tbx_tick_t, *p_tbx_tick_t;


#if defined(CRLIBM_TYPECPU_ITANIUM) && defined(__GNUC__) && !defined(__INTEL_COMPILER)
 
#define TBX_GET_TICK(t) \
   __asm__ __volatile__("mov %0=ar.itc" : "=r"((t).tick) :: "memory")
#define TBX_TICK_RAW_DIFF(t1, t2)    ((t2).tick - (t1).tick)


#elif defined(CRLIBM_TYPECPU_ITANIUM) && defined(__HPCC__)
/* Never tested, currently nobody defines __HPCC__.
  It should work on HPUX machines using HP compiler */
#include<machine/sys/inline.h>/* to move up of course */
#define TBX_GET_TICK(t) \
   t.tick=_Asm_mov_from_ar(_AREG_ITC)
#define TBX_TICK_RAW_DIFF(t1, t2)    ((t2).tick - (t1).tick )


#elif defined(CRLIBM_TYPECPU_ITANIUM) && defined(__INTEL_COMPILER)
#define TBX_GET_TICK(t)    t.tick=__getReg(_IA64_REG_AR_ITC)   
#define TBX_TICK_RAW_DIFF(t1, t2)    ((t2).tick - (t1).tick)



#elif (defined(CRLIBM_TYPECPU_AMD64) || defined(CRLIBM_TYPECPU_X86)) \
    && defined(__GNUC__) && !defined(__INTEL_COMPILER)
#define TBX_GET_TICK(time)                \
        __asm__ __volatile__(             \
                "xorl %%eax,%%eax\n\t"    \
                "cpuid\n\t"               \
                "rdtsc\n\t"               \
                "movl %%eax,(%0)\n\t"     \
                "movl %%edx,4(%0)\n\t"    \
                "xorl %%eax,%%eax\n\t"    \
                "cpuid\n\t"               \
                : /* no output */         \
                : "S"(&time.tick)         \
                : "eax", "ebx", "ecx", "edx", "memory")
#define TBX_TICK_RAW_DIFF(t1, t2)    ((t2).tick - (t1).tick)


#elif defined(CRLIBM_TYPECPU_ALPHA)  && defined(__GNUC__)
#define TBX_GET_TICK(t) \
   __asm__ volatile("rpcc %0\n\t" : "=r"((t).tick))
#define TBX_TICK_RAW_DIFF(t1, t2) \
   (((t2).tick & 0xFFFFFFFF) - ((t1).tick & 0xFFFFFFFF))


#elif defined(CRLIBM_TYPECPU_SPARC)  && defined(__GNUC__)
#define TBX_GET_TICK(t) \
    (t).tick = gethrtime()
#define TBX_TICK_RAW_DIFF(t1, t2) \
   ((t2).tick  - (t1).tick)



#else  /* generic case */

#define TIMING_USES_GETTIMEOFDAY 1  /* inform the test functions that timings are inaccurate */
#define TBX_GET_TICK(t) \
   gettimeofday(&(t).timev, 0)
#define TBX_TICK_RAW_DIFF(t1, t2) \
   ((t2.timev.tv_sec * 1000000L + t2.timev.tv_usec) - \
    (t1.timev.tv_sec * 1000000L + t1.timev.tv_usec))

#endif

#define TBX_TICK_DIFF(t1, t2) (TBX_TICK_RAW_DIFF(t1, t2) - tbx_residual + 1)
#define TBX_TIMING_DELAY(t1, t2) tbx_tick2usec(TBX_TICK_DIFF(t1, t2))

extern unsigned long long tbx_residual;
extern tbx_tick_t         tbx_new_event;
extern tbx_tick_t         tbx_last_event;

#endif /* TBX_TIMING_H */
