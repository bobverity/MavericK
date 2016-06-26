//
//  MavericK
//  MCMC_TI_noAdmixture.cpp
//
//  Created: Bob on 22/06/2016
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "MCMC_TI_noAdmixture.h"

using namespace std;

//------------------------------------------------
// MCMC_TI_noAdmixture::
// constructor for class containing all elements required for MCMC under no-admixture model
MCMC_TI_noAdmixture::MCMC_TI_noAdmixture(globals &globals, int _Kindex, int _rungs) {
    
    // copy some values over from globals object
    Kindex = _Kindex;
    K = globals.Kmin+Kindex;
    burnin = globals.thermodynamicBurnin;
    samples = globals.thermodynamicSamples;
    thinning = globals.thermodynamicThinning;
    
    // create a chain for each rung
    rungs = _rungs;
    betaVec = vector<double>(rungs);
    chainVec = vector<chain_noAdmixture>(rungs);
    for (int rung=0; rung<rungs; rung++) {
        betaVec[rung] = rung/double(rungs-1);
        chainVec[rung] = chain_noAdmixture(globals, K, betaVec[rung]);
    }
    
    // likelihoods for each rung
    logLikeGroup_sum = vector<double>(rungs);
    logLikeGroup_sumSquared = vector<double>(rungs);
    logLikeGroup_store = vector< vector<double> >(rungs, vector<double>(samples));
    
    // autocorrelation for each rung
    autoCorr = vector<double>(rungs);
    ESS = vector<double>(rungs);
    
    // TI point estimates for each rung
    TIpoint_mean = vector<double>(rungs);
    TIpoint_var = vector<double>(rungs);
    TIpoint_SE = vector<double>(rungs);
    
    // overall TI point estimate and SE
    logEvidence_TI = 0;
    logEvidence_TI_var = 0;
    logEvidence_TI_SE = 0;
    
}

//------------------------------------------------
// MCMC_TI_noAdmixture::
// perform complete MCMC under no-admixture model
void MCMC_TI_noAdmixture::perform_MCMC() {
    
    // reset chains
    for (int rung=0; rung<rungs; rung++) {
        chainVec[rung].reset(true);
    }
    
    // perform MCMC
    int thinSwitch = 1;
    for (int rep=0; rep<(burnin+samples); rep++) {
        
        // thinning loop (becomes active after burn-in)
        for (int thin=0; thin<thinSwitch; thin++) {
            
            // update group allocation of all individuals in all rungs
            for (int rung=0; rung<rungs; rung++) {
                chainVec[rung].group_update();
            }
            
        }
        if (rep==burnin)
            thinSwitch = thinning;
        
        // loop over rungs
        for (int rung=0; rung<rungs; rung++) {
        
            // calculate marginal likelihood
            chainVec[rung].d_logLikeGroup();
            
            // add likelihoods to running sums
            if (rep>=burnin) {
                logLikeGroup_sum[rung] += chainVec[rung].logLikeGroup;
                logLikeGroup_sumSquared[rung] += chainVec[rung].logLikeGroup*chainVec[rung].logLikeGroup;
                logLikeGroup_store[rung][rep-burnin] = chainVec[rung].logLikeGroup;
            }
            
        }
        
    } // end of MCMC
    
    // process likelihoods for each rung
    for (int rung=0; rung<rungs; rung++) {
        
        // calculate autocorrelation
        autoCorr[rung] = calculateAutoCorr(logLikeGroup_store[rung]);
        ESS[rung] = samples/autoCorr[rung];
        
        // calculate mean and SE of log-likelihood
        TIpoint_mean[rung] = logLikeGroup_sum[rung]/double(samples);
        TIpoint_var[rung] = logLikeGroup_sumSquared[rung]/double(samples) - TIpoint_mean[rung]*TIpoint_mean[rung];
        if (TIpoint_var[rung]<0) {   // avoid arithmetic underflow
            TIpoint_var[rung] = 0;
        }
        TIpoint_SE[rung] = sqrt(TIpoint_var[rung]/ESS[rung]);
        
    }
    
    // calculate overall TI point estimate and SE. Note that simply summing the variance by looping over rungs will give the wrong answer, as it fails to account for the fact that variance(2*A) is not 2*variance(A) but rather 4*variance(A).
    logEvidence_TI = 0.5*TIpoint_mean[0]*(betaVec[1]-betaVec[0]) + 0.5*TIpoint_mean[rungs-1]*(betaVec[rungs-1]-betaVec[rungs-2]);
    logEvidence_TI_var = 0.25*TIpoint_var[0]/ESS[0]*pow(betaVec[1]-betaVec[0],2) + 0.25*TIpoint_var[rungs-1]/ESS[rungs-1]*pow(betaVec[rungs-1]-betaVec[rungs-2],2);
    if (rungs>2) {
        for (int rung=1; rung<(rungs-1); rung++) {
            logEvidence_TI += TIpoint_mean[rung]*(betaVec[rung] - betaVec[rung-1]);
            logEvidence_TI_var += TIpoint_var[rung]/ESS[rung]*pow(betaVec[rung]-betaVec[rung-1],2);
        }
    }
    logEvidence_TI_SE = sqrt(logEvidence_TI_var);
    
}

