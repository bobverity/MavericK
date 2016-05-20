//
//  MavericK
//  EM_algorithm.h
//
//  Created: Bob on 25/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  These functions implement the expectation-maximization (EM) algorithm to find maximum likelihood allele frequencies and admixture proportions under both the with- and wothout-admixture models. The final maximum likelihood values are used by various model comparison statistics.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__EM_algorithm__
#define __Maverick1_0__EM_algorithm__

#include <iostream>

#include "globals.h"
#include "probability.h"
#include "misc.h"

//------------------------------------------------
// EM algorithm under no-admixture model
void EM_noAdmix(globals &globals, int Kindex);

//------------------------------------------------
// EM algorithm under admixture model
void EM_admix(globals &globals, int Kindex);

#endif
