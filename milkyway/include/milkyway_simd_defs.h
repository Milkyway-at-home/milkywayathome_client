
#ifndef _MILKYWAY_SIMD_DEFS_H_
#define _MILKYWAY_SIMD_DEFS_H_

#define DZERO           (mw_set1_pd(0.0))
#define DONEHALF        (mw_set1_pd(0.5))
#define DONEHALF_NEG    (mw_set1_pd(-0.5))
#define DONE            (mw_set1_pd(1.0))
#define DONE_NEG        (mw_set1_pd(-1.0))
#define DTWO            (mw_set1_pd(2.0))
#define DTHREE          (mw_set1_pd(3.0))
#define DFOUR           (mw_set1_pd(4.0))
#define DFIVE           (mw_set1_pd(5.0))
#define LOG2ME          (mw_set1_pd(1.4426950408889634073599))
#define L2U             0.69314718055966295651160180568695068359375
#define L2L             0.28235290563031577122588448175013436025525412068e-12

#endif /* _MILKYWAY_SIMD_DEFS_H_ */

