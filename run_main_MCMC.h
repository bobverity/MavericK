//
//  MavericK
//  run_main_MCMC.h
//
//  Created: Bob on 23/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  These functions implement the main MCMC, in which allocations of all individuals or gene copies are drawn from the posterior distribution. Note that the main MCMC does not implement thermodynamic integration, but rather is used to obtain Q-matrix output and to calculate statistics such as the "Structure" estimator L_K.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__run_main_MCMC__
#define __Maverick1_0__run_main_MCMC__

#include <iostream>
#include "globals.h"
#include "probability.h"
#include "misc.h"
#include "MCMCobject_admixture.h"
#include "MCMC_main_noAdmixture.h"

//------------------------------------------------
// main MCMC under no-admixture model, repeated multiple times
void run_main_MCMC_noAdmixture(globals &globals, int Kindex);

//------------------------------------------------
// main MCMC under admixture model, repeated multiple times
void mainMCMC_admixture(globals &globals, int Kindex);


#endif
