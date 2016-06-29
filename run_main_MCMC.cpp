//
//  MavericK
//  run_main_MCMC.cpp
//
//  Created: Bob on 23/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "run_main_MCMC.h"

using namespace std;

//------------------------------------------------
// main Structure MCMC under admixture model, repeated multiple times
void mainMCMC_admixture(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // define MCMC object
    MCMCobject_admixture mainMCMC(globals, Kindex, globals.mainBurnin, globals.mainSamples, 1, 1.0);
        
    // perform MCMC
    mainMCMC.reset(true);
    mainMCMC.perform_MCMC(globals, true, false, true, globals.outputLikelihood_on, globals.outputPosteriorGrouping_on, 1);
    
    // save Qmatrix values
    for (int i=0; i<globals.geneCopies; i++) {
        for (int k=0; k<K; k++) {
            globals.Qmatrix_gene[Kindex][i][k] = mainMCMC.Qmatrix_gene[i][k];
        }
    }
    for (int i=0; i<globals.n; i++) {
        for (int k=0; k<K; k++) {
            globals.Qmatrix_ind[Kindex][i][k] = mainMCMC.Qmatrix_ind[i][k];
        }
    }
    if (globals.outputQmatrix_pop_on) {
        for (int i=0; i<int(globals.uniquePops.size()); i++) {
            for (int k=0; k<K; k++) {
                globals.Qmatrix_pop[Kindex][i][k] = mainMCMC.Qmatrix_pop[i][k];
            }
        }
    }
    
    // save harmonic mean of marginal likelihoods
    globals.logEvidence_harmonic[Kindex] = mainMCMC.harmonic;
    
    // calculate Structure estimator from joint likelihoods
    globals.structure_loglike_mean[Kindex] = mainMCMC.logLikeJoint_sum/double(globals.mainSamples);
    globals.structure_loglike_var[Kindex] = mainMCMC.logLikeJoint_sumSquared/double(globals.mainSamples) - globals.structure_loglike_mean[Kindex]*globals.structure_loglike_mean[Kindex];
    globals.logEvidence_structure[Kindex] = globals.structure_loglike_mean[Kindex] - 0.5*globals.structure_loglike_var[Kindex];
    
}