#include <stdio.h>
#include <stdlib.h>
#include "evaluator.h"

void init_simple_evaluator(double (*lf)(double*)) {
	evaluate = lf;
}
