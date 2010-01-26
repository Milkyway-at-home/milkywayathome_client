//doesn't appear to work as a macro
#define SQRT							\
  //manually convert double to float (ATI Hardware can't do it)		\
  unsigned long rgl = as_ulong(rg);					\
  int sign, exponent;							\
  float value;								\
  sign = (int) ((rgl & 0x8000000000000000) >> 63);			\
  exponent = (int) ((rgl & 0x7FE0000000000000) >> 52);			\
  value = (rgl & 0x001FFFFFFFFFFFFF) + 0x0010000000000000;		\
  value = ldexp(value, exponent - 1023 - 52);				\
  if (sign)								\
    value = -value;							\
  //rsqrt(y) -> http://en.wikipedia.org/wiki/Fast_inverse_square_root	\
  double x  = (double) rsqrt(value);					\
  //fsqrtd								\
  x = x * (3.0 - (rg)*(x*x));						\
  double res = x * (rg);						\
  rg = res * (0.75 - 0.0625*(res*x));					\
  
