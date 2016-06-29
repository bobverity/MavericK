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
    n = globals.n;
    uniquePops = globals.uniquePops;
    outputQmatrix_pop_on = globals.outputQmatrix_pop_on;
    burnin = globals.mainBurnin;
    samples = globals.mainSamples;
    
    // create a chain for each rung
    rungs = _rungs;
    betaVec = vector<double>(rungs);
    chainVec = vector<chain_noAdmixture>(rungs);
    for (int rung=0; rung<rungs; rung++) {
        betaVec[rung] = rung/double(rungs-1);
        chainVec[rung] = chain_noAdmixture(globals, K, betaVec[rung]);
    }
    acceptanceRate = vector<double>(rungs-1);
    
    // likelihoods for each rung
    logLikeGroup_sum = vector<double>(rungs);
    logLikeGroup_sumSquared = vector<double>(rungs);
    logLikeGroup_store = vector< vector<double> >(rungs, vector<double>(samples));
    logLikeJoint_sum = 0;
    logLikeJoint_sumSquared = 0;
    harmonic = log(double(0));
    
    // initialise Qmatrices
    Qmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_pop = vector< vector<double> >(uniquePops.size(), vector<double>(K));
    
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
void MCMC_TI_noAdmixture::perform_MCMC(globals &globals, bool drawAlleleFreqs, bool fixLabels, bool outputLikelihood, bool outputPosteriorGrouping, int mainRep) {
    
    // reset chains
    for (int rung=0; rung<rungs; rung++) {
        chainVec[rung].reset(true);
    }
    
    acceptanceRate = vector<double>(rungs-1);
    
    // perform MCMC
    for (int rep=0; rep<(burnin+samples); rep++) {
            
        // update group allocation of all individuals in all rungs
        for (int rung=0; rung<rungs; rung++) {
            chainVec[rung].group_update();
        }
        
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
        
        // Metropolis-coupling
        // loop over rungs, starting with the second hottest chain and moving to the cold chain
        for (int rung=1; rung<rungs; rung++) {
            
            // get log-likelihoods of hot and cold chain in the comparison
            double coldLikelihood = chainVec[rung].logLikeGroup;
            double hotLikelihood = chainVec[rung-1].logLikeGroup;
            
            // calculate acceptance ratio (still in log space)
            double acceptance = (hotLikelihood*betaVec[rung] + coldLikelihood*betaVec[rung-1]) - (coldLikelihood*betaVec[rung] + hotLikelihood*betaVec[rung-1]);
            
            // accept or reject move
            double rand1 = runif1();
            if (log(rand1)<acceptance && 1==1) {
                
                spareChain = chainVec[rung];
                Qmatrix_spare = chainVec[rung].logQmatrix_ind;
                Qmatrix_spare2 = chainVec[rung].logQmatrix_ind_running;
                
                chainVec[rung] = chainVec[rung-1];
                chainVec[rung-1] = spareChain;
                
                chainVec[rung].beta = betaVec[rung];
                chainVec[rung-1].beta = betaVec[rung-1];
                
                chainVec[rung].logQmatrix_ind = Qmatrix_spare;
                chainVec[rung].logQmatrix_ind_running = Qmatrix_spare2;
                
                acceptanceRate[rung-1] += 1.0/double(burnin+samples);
                
            }
            
        }
        
        // print to junk file
        if (rep>=burnin) {
            for (int rung=0; rung<rungs; rung++) {
                globals.junk_fileStream << chainVec[rung].logLikeGroup << "\t";
            }
            globals.junk_fileStream << "\n";
        }
        
        // if fix label-switching problem
        if (fixLabels) {
            
            // fix label-switching problem
            chainVec[rungs-1].chooseBestLabelPermutation(globals);
            
            // add logQmatrix_ind_new to logQmatrix_ind_running
            chainVec[rungs-1].updateQmatrix();
            
            // store Qmatrix values if no longer in burn-in phase
            if (rep>=burnin)
                chainVec[rungs-1].storeQmatrix();
        }
        
        // optionally draw allele frequencies and calculate joint likelihood
        if (drawAlleleFreqs) {
            chainVec[rungs-1].drawFreqs();
            chainVec[rungs-1].d_logLikeJoint();
        }
        
        // add likelihoods to running sums
        if (rep>=burnin) {
            
            for (int rung=0; rung<rungs; rung++) {
                logLikeGroup_sum[rung] += chainVec[rung].logLikeGroup;
                logLikeGroup_sumSquared[rung] += chainVec[rung].logLikeGroup*chainVec[rung].logLikeGroup;
            }
            harmonic = logSum(harmonic, -chainVec[rungs-1].logLikeGroup);
            
            if (drawAlleleFreqs) {
                logLikeJoint_sum += chainVec[rungs-1].logLikeJoint;
                logLikeJoint_sumSquared += chainVec[rungs-1].logLikeJoint*chainVec[rungs-1].logLikeJoint;
            }
        }
        
        // write to outputLikelihoods file
        if (outputLikelihood) {
            globals.outputLikelihood_fileStream << K << "," << mainRep+1 << "," << rep-burnin+1 << "," << chainVec[rungs-1].logLikeGroup << "," << chainVec[rungs-1].logLikeJoint << "\n";
            globals.outputLikelihood_fileStream.flush();
        }
        
        // write to outputPosteriorGrouping file
        if (outputPosteriorGrouping) {
            globals.outputPosteriorGrouping_fileStream << K << "," << mainRep+1 << "," << rep-burnin+1;
            for (int i=0; i<n; i++) {
                globals.outputPosteriorGrouping_fileStream << "," << chainVec[rungs-1].group[i];
            }
            globals.outputPosteriorGrouping_fileStream << "\n";
            globals.outputPosteriorGrouping_fileStream.flush();
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

