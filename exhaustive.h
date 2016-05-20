//
//  MavericK
//  exhaustive.h
//
//  Created: Bob on 23/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  These functions search all possible partitions of the data to compute the exact model evidence. This is only possible for extremely small data sets, on the order of n=10 individuals under the without-admixture model, or even fewer under the admixture model.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__exhaustive__
#define __Maverick1_0__exhaustive__

#include <iostream>
#include "globals.h"
#include "misc.h"

//------------------------------------------------
// exhaustive analysis under no-admixture model
void exhaustive_noAdmix(globals &globals, int Kindex);

//------------------------------------------------
// exhaustive analysis under admixture model (alpha fixed or variable)
void exhaustive_admix(globals &globals, int Kindex);

//------------------------------------------------
// exhaustive analysis under admixture model for given alpha
double exhaustive_admix_fixedAlpha(globals &globals, int Kindex, double alpha);

#endif
