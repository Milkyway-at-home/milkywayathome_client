#ifndef ST_CNUM_H
#define ST_CNUM_H

typedef struct {
	double real;
	double imagine;
} cnum;

/* complex print function */
void PrintCnum(cnum *x, int n);

/* complex add functions */
cnum CnumAdd(cnum c1, cnum c2);
cnum CnumAddD(cnum c1, double d);

/* complex subtraction functions */
cnum CnumSub(cnum c1, cnum c2);
cnum CnumSubD(cnum c1, double d);
cnum CnumDSub(double d, cnum c1);

/* complex multiplication functions */
cnum CnumMult(cnum c1, cnum c2);
cnum CnumMultD(cnum c1, double d);

/* complex division functions */
cnum CnumDiv(cnum c1, cnum c2);
cnum CnumDivD(cnum c1, double d);
cnum CnumDDiv(double x, cnum c1);

/* complex square root and cube root functions */
cnum CnumSqrt(cnum c1);
cnum CnumCbrt(cnum c1, int verb);

#endif
