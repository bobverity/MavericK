//
//  MavericK
//  MCMC_main_noAdmixture.cpp
//
//  Created: Bob on 22/06/2016
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "MCMC_main_noAdmixture.h"

using namespace std;

//------------------------------------------------
// MCMC_main_noAdmixture::
// constructor for class containing all elements required for MCMC under no-admixture model
MCMC_main_noAdmixture::MCMC_main_noAdmixture(globals &globals, int _Kindex, int _rungs) {
    
    // copy some values over from globals object
    Kindex = _Kindex;
    K = globals.Kmin+Kindex;
    n = globals.n;
    uniquePops = globals.uniquePops;
    outputQmatrix_pop_on = globals.outputQmatrix_pop_on;
    burnin = globals.mainBurnin;
    samples = globals.mainSamples;
    thinning = globals.mainThinning;
    
    // create a chain for each rung used in Metropolis-coupling
    //rungs = _rungs;
    rungs = 20;
    betaVec = vector<double>(rungs);
    chainVec = vector<chain_noAdmixture>(rungs);
    for (int rung=0; rung<rungs; rung++) {
        betaVec[rung] = rung/double(rungs-1);
        chainVec[rung] = chain_noAdmixture(globals, K, betaVec[rung]);
    }
    acceptanceRate = vector<double>(rungs-1);
    
    // likelihoods
    logLikeGroup_sum = 0;
    logLikeGroup_sumSquared = 0;
    logLikeJoint_sum = 0;
    logLikeJoint_sumSquared = 0;
    harmonic = log(double(0));
    
    // initialise Qmatrices
    Qmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_pop = vector< vector<double> >(uniquePops.size(), vector<double>(K));
    
}

//------------------------------------------------
// MCMC_main_noAdmixture::
// reset all chains, with option for whether to reset running estimate of Q-matrix
void MCMC_main_noAdmixture::reset(bool reset_Qmatrix_running) {
    
    // reset chains
    for (int rung=0; rung<rungs; rung++) {
        chainVec[rung].reset(reset_Qmatrix_running);
    }
    
}

//------------------------------------------------
// MCMC_main_noAdmixture::
// perform complete MCMC under no-admixture model
void MCMC_main_noAdmixture::perform_MCMC(globals &globals, bool drawAlleleFreqs, bool fixLabels, bool outputLikelihood, bool outputPosteriorGrouping, int mainRep) {
    
    acceptanceRate = vector<double>(rungs-1);
    
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
        
        // calculate marginal likelihood in all rungs
        for (int rung=0; rung<rungs; rung++) {
            chainVec[rung].d_logLikeGroup();
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
        /*
        for (int rung=0; rung<rungs; rung++) {
            
            // fix label-switching problem
            if (rung==(rungs-1))
                chainVec[rung].chooseBestLabelPermutation(globals);
            
            // add logQmatrix_ind_new to logQmatrix_ind_running
            if (rung==(rungs-1))
                chainVec[rung].updateQmatrix();
            
            // store Qmatrix values if no longer in burn-in phase
            if (rep>=burnin && rung==(rungs-1))
                chainVec[rung].storeQmatrix();
        }
        */
        
        // optionally draw allele frequencies and calculate joint likelihood
        if (drawAlleleFreqs) {
            chainVec[rungs-1].drawFreqs();
            chainVec[rungs-1].d_logLikeJoint();
        }
        
        // add likelihoods to running sums
        if (rep>=burnin) {
            logLikeGroup_sum += chainVec[rungs-1].logLikeGroup;
            logLikeGroup_sumSquared += chainVec[rungs-1].logLikeGroup*chainVec[rungs-1].logLikeGroup;
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

    
    // finish off Qmatrices
    if (fixLabels) {
        
        // finish off individual level Qmatrix
        for (int i=0; i<n; i++) {
            for (int k=0; k<K; k++) {
                Qmatrix_ind[i][k] = exp(chainVec[rungs-1].logQmatrix_ind[i][k] - log(double(samples)));
            }
        }
        
        // calculate population level Qmatrices
        if (outputQmatrix_pop_on) {
            for (int i=0; i<n; i++) {
                for (int k=0; k<K; k++) {
                    Qmatrix_pop[globals.pop_index[i]][k] += Qmatrix_ind[i][k];
                }
            }
            for (int i=0; i<int(uniquePops.size()); i++) {
                for (int k=0; k<K; k++) {
                    Qmatrix_pop[i][k] /= double(globals.uniquePop_counts[i]);
                }
            }
        }
    } // end if fixLabels
    
    // finish off harmonic mean
    harmonic = log(double(samples))-harmonic;
    
    //printVector(acceptanceRate);
}

