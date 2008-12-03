#include <stdio.h>
#include <stdlib.h>
#include "evaluator.h"

double (*evaluate)(double*);

void init_simple_evaluator(double (*lf)(double*)) {
	evaluate = lf;
}
