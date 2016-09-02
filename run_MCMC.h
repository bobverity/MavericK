//
//  MavericK
//  run_MCMC.h
//
//  Created: Bob on 23/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Carries out thermodynamic integration for models both with- and without-admixture.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__run_MCMC__
#define __Maverick1_0__run_MCMC__

#include <iostream>
#include "globals.h"
#include "MCMCobject_admixture.h"
#include "misc.h"
#include "MCMC_noAdmixture.h"

//------------------------------------------------
// thermodynamic integral estimator for no-admixture model
void run_MCMC_noAdmixture(globals &globals, int Kindex);

//------------------------------------------------
// thermodynamic integral estimator for admixture model
void TI_admixture(globals &globals, int Kindex);

#endif
