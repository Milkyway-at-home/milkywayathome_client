/*
 *  parameters.h
 *  Astronomy
 *
 *  Created by Travis Desell on 2/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef BOINC_ASTRONOMY_PARAMETERS_H
#define BOINC_ASTRONOMY_PARAMETERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parameters.h"

int	boinc_read_astronomy_parameters(const char* file, ASTRONOMY_PARAMETERS *ap);
int	boinc_write_astronomy_parameters(const char* file, ASTRONOMY_PARAMETERS *ap);

#endif
