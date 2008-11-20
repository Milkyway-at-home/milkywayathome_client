#ifndef ASTRONOMY_EVALUATION_H
#define ASTRONOMY_EVALUATION_H

#include "evaluation_state.h"
#include "parameters.h"
#include "star_points.h"

int	calculate_integrals(ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, STAR_POINTS* sp);
int	calculate_likelihood(ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, STAR_POINTS* sp);

#endif


