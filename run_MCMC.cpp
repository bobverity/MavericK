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
    mainMCMC.perform_MCMC(globals, globals.outputLikelihood_on, globals.outputPosteriorGrouping_on);
    
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
    
    // calculate Structure estimator from joint likelihoods
    globals.structure_loglike_mean[Kindex] = mainMCMC.logLikeJoint_sum/double(globals.mainSamples);
    globals.structure_loglike_var[Kindex] = mainMCMC.logLikeJoint_sumSquared/double(globals.mainSamples) - globals.structure_loglike_mean[Kindex]*globals.structure_loglike_mean[Kindex];
    globals.logEvidence_structure[Kindex] = globals.structure_loglike_mean[Kindex] - 0.5*globals.structure_loglike_var[Kindex];
    
    // save results to global variable
    globals.TIpoint_mean[Kindex] = mainMCMC.TIpoint_mean;
    globals.TIpoint_var[Kindex] = mainMCMC.TIpoint_var;
    globals.TIpoint_SE[Kindex] = mainMCMC.TIpoint_SE;
    globals.logEvidence_TI[Kindex] = mainMCMC.logEvidence_TI;
    globals.logEvidence_TI_SE[Kindex] = mainMCMC.logEvidence_TI_SE;
    
}

//------------------------------------------------
// run main MCMC for admixture model
void TI_admixture(globals &globals, int Kindex) {
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
    
    // set up beta vector
    vector<double> betaVec(globals.mainRungs);
    for (int TIrep=0; TIrep<globals.mainRungs; TIrep++)
        betaVec[TIrep] = double(TIrep)/(globals.mainRungs-1);
    
    // carry out thermodynamic integration
    for (int TIrep=0; TIrep<globals.mainRungs; TIrep++) {
        double beta = betaVec[TIrep];
        
        char * buffer = new char[8];
        sprintf(buffer, "%.2f", beta);
        string s = buffer;
        coutAndLog("  power = "+s+"\n", globals.outputLog_on, globals.outputLog_fileStream);
        
        // define MCMC object
        MCMCobject_admixture TI_MCMC(globals, Kindex, globals.mainBurnin, globals.mainSamples, 1, beta);
        
        // perform MCMC
        TI_MCMC.reset(true);
		TI_MCMC.perform_MCMC(globals, false, true, false, false, false, TIrep);
        
        // calculate autocorrelation
        double autoCorr = calculateAutoCorr(TI_MCMC.logLikeGroup_store);
        double ESS = globals.mainSamples/autoCorr;
        
        // save results of this power
        globals.TIpoint_mean[Kindex][TIrep] = TI_MCMC.logLikeGroup_sum/double(globals.mainSamples);
        globals.TIpoint_var[Kindex][TIrep] = TI_MCMC.logLikeGroup_sumSquared/double(globals.mainSamples) - globals.TIpoint_mean[Kindex][TIrep]*globals.TIpoint_mean[Kindex][TIrep];
        globals.TIpoint_SE[Kindex][TIrep] = sqrt(globals.TIpoint_var[Kindex][TIrep]/ESS);
        
    }
    
    // calculate thermodynamic integral estimate. Note that at this stage we assume equally spaced rungs.
    double Dsum = 0.5*globals.TIpoint_mean[Kindex][0] + 0.5*globals.TIpoint_mean[Kindex][globals.mainRungs-1];
    double Vsum = 0.25*globals.TIpoint_SE[Kindex][0]*globals.TIpoint_SE[Kindex][0] + 0.25*globals.TIpoint_SE[Kindex][globals.mainRungs-1]*globals.TIpoint_SE[Kindex][globals.mainRungs-1];
    if (globals.mainRungs>2) {
        for (int TIrep=1; TIrep<(globals.mainRungs-1); TIrep++) {
            Dsum += globals.TIpoint_mean[Kindex][TIrep];
            Vsum += globals.TIpoint_SE[Kindex][TIrep]*globals.TIpoint_SE[Kindex][TIrep];
        }
    }
    Dsum /= double(globals.mainRungs-1);
    Vsum /= double(globals.mainRungs-1)*double(globals.mainRungs-1);
    globals.logEvidence_TI[Kindex] = Dsum;
    globals.logEvidence_TI_SE[Kindex] = sqrt(Vsum);
    
}