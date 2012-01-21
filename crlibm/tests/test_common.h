extern double zero;

void test_init(/* pointers to returned value */
	       double (**randfun_perf)(), 
	       double (**randfun_soaktest)(), 
	       double (**testfun_crlibm)(), 
	       int    (**testfun_mpfr)  (),
	       double (**testfun_libultim)   (),
	       double (**testfun_libmcr)  (),
	       double (**testfun_libm)  (),
	       double* worst_case,
	       /* arguments */
	       char *func_name,
	       char *rnd_mode);

int rand_int(void);

double rand_generic(void);

double rand_double(void);

double rand_for_exp(void);

double rand_for_exp_perf(void);

double rand_for_log(void);

double rand_for_exp_normal(void);

double rand_for_expm1_soaktest(void);

