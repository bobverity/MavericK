//
//  MavericK
//  run_MCMC.cpp
//
//  Created: Bob on 23/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "run_MCMC.h"

using namespace std;

//------------------------------------------------
// run main MCMC for no-admixture model
void run_MCMC_noAdmixture(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // special case if K==1
    if (K==1) {
        for (int TIrep=0; TIrep<globals.mainRungs; TIrep++) {
            globals.TIpoint_mean[Kindex][TIrep] = globals.logEvidence_exhaustive[Kindex];
            globals.TIpoint_var[Kindex][TIrep] = 0;
            globals.TIpoint_SE[Kindex][TIrep] = 0;
        }
        globals.logEvidence_TI[Kindex] = globals.logEvidence_exhaustive[Kindex];
        globals.logEvidence_TI_SE[Kindex] = 0;
        return;
    }
    
    // perform MCMC
    MCMC_noAdmixture mainMCMC(globals, Kindex, globals.mainRungs);
    mainMCMC.perform_MCMC(globals);
    
    // save Qmatrix values
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
    
    // save Structure estimator
    globals.structure_loglike_mean[Kindex] = mainMCMC.structure_loglike_mean;
    globals.structure_loglike_var[Kindex] = mainMCMC.structure_loglike_var;
    globals.logEvidence_structure[Kindex] = mainMCMC.logEvidence_structure;
    
    // save TI results to global variable
    globals.TIpoint_mean[Kindex] = mainMCMC.TIpoint_mean;
    globals.TIpoint_var[Kindex] = mainMCMC.TIpoint_var;
    globals.TIpoint_SE[Kindex] = mainMCMC.TIpoint_SE;
    globals.logEvidence_TI[Kindex] = mainMCMC.logEvidence_TI;
    globals.logEvidence_TI_SE[Kindex] = mainMCMC.logEvidence_TI_SE;
}

//------------------------------------------------
// run main MCMC for admixture model
void run_MCMC_admixture(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // special case if K==1
    if (K==1) {
        for (int TIrep=0; TIrep<globals.mainRungs; TIrep++) {
            globals.TIpoint_mean[Kindex][TIrep] = globals.logEvidence_exhaustive[Kindex];
            globals.TIpoint_var[Kindex][TIrep] = 0;
            globals.TIpoint_SE[Kindex][TIrep] = 0;
        }
        globals.logEvidence_TI[Kindex] = globals.logEvidence_exhaustive[Kindex];
        globals.logEvidence_TI_SE[Kindex] = 0;
        return;
    }
    
    // perform MCMC
    MCMC_admixture mainMCMC(globals, Kindex, globals.mainRungs);
    mainMCMC.perform_MCMC(globals);
    
    // save Qmatrix values
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
    
    // save Structure estimator
    globals.structure_loglike_mean[Kindex] = mainMCMC.structure_loglike_mean;
    globals.structure_loglike_var[Kindex] = mainMCMC.structure_loglike_var;
    globals.logEvidence_structure[Kindex] = mainMCMC.logEvidence_structure;
    
    // save TI results to global variable
    globals.TIpoint_mean[Kindex] = mainMCMC.TIpoint_mean;
    globals.TIpoint_var[Kindex] = mainMCMC.TIpoint_var;
    globals.TIpoint_SE[Kindex] = mainMCMC.TIpoint_SE;
    globals.logEvidence_TI[Kindex] = mainMCMC.logEvidence_TI;
    globals.logEvidence_TI_SE[Kindex] = mainMCMC.logEvidence_TI_SE;
    
}
