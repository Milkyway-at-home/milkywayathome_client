
#ifndef _MILKYWAY_SIMD_H_
#define _MILKYWAY_SIMD_H_

#if defined(__AVX__) && 0
  /* Needs work */
  #include "milkyway_avx_intrin.h"
#elif defined(__SSE2__) || defined(__SSE3__) || defined(__SSSE3__) || defined(__SSE41__)
  #include "milkyway_sse2_intrin.h"
#endif


#define ROUND_DOWN(x, s) ((x) & ~((s)-1))
union d2
{
    int i[4];
    unsigned int u[4];
    long int lu[2];
    double d[2];
    __m128d m;
    __m128i mi;
};

#endif /* _MILKYWAY_SIMD_H_ */

