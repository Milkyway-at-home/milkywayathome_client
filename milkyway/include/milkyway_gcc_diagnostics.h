
#ifndef _MILKYWAY_GCC_DIAGNOSTICS_H_
#define _MILKYWAY_GCC_DIAGNOSTICS_H_

/*
  https://svn.boost.org/trac/boost/wiki/Guidelines/WarningsGuidelines
*/

#ifndef __INTEL_COMPILER
  #if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 402
    #define GCC_DIAG_STR(s) #s
    #define GCC_DIAG_JOINSTR(x,y) GCC_DIAG_STR(x ## y)
    #define GCC_DIAG_DO_PRAGMA(x) _Pragma (#x)
    #define GCC_DIAG_PRAGMA(x) GCC_DIAG_DO_PRAGMA(GCC diagnostic x)
    #if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
      #define GCC_DIAG_OFF(x) GCC_DIAG_PRAGMA(push) \
          GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x))
      #define GCC_DIAG_ON(x) GCC_DIAG_PRAGMA(pop)
    #else
      #define GCC_DIAG_OFF(x) GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x))
      #define GCC_DIAG_ON(x)  GCC_DIAG_PRAGMA(warning GCC_DIAG_JOINSTR(-W,x))
    #endif
  #else
    #define GCC_DIAG_OFF(x)
    #define GCC_DIAG_ON(x)
  #endif
#else
  # define GCC_DIAG_OFF(x)
  # define GCC_DIAG_ON(x)
#endif /*__INTEL_COMPILER */

#endif /* _MILKYWAY_GCC_DIAGNOSTICS_H_ */

