/* Vendored CORE-MATH correctly-rounded functions (binary64) that crlibm
 * does not provide. https://core-math.gitlabpages.inria.fr/  (MIT-style).
 * Built with the normal milkyway flags (no -mfma); the internal fma()
 * resolves to the libm software fma so no hardware FMA instruction is
 * emitted (runs on devices without FMA). */
#ifndef _MW_COREMATH_H_
#define _MW_COREMATH_H_
#ifdef __cplusplus
extern "C" {
#endif
double cr_cbrt(double x);
double cr_erf(double x);
#ifdef __cplusplus
}
#endif
#endif
