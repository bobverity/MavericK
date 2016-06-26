//
//  MavericK
//  MCMC_main_noAdmixture.h
//
//  Created: Bob on 22/06/2016
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Defines a class that can be used to carry out MCMC under the without-admixture model.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__MCMC_main_noAdmixture__
#define __Maverick1_0__MCMC_main_noAdmixture__

#include <iostream>
#include "globals.h"
#include "probability.h"
#include "misc.h"
#include "Hungarian.h"
#include "chain_noAdmixture.h"

//------------------------------------------------
// class containing all elements required for MCMC under no-admixture model
class MCMC_main_noAdmixture {
    
public:
    
    // PUBLIC OBJECTS
    
    //chain_noAdmixture chain;
    
    // basic quantities
    int Kindex;
    int K;
    int n;
    std::vector<std::string> uniquePops;
    bool outputQmatrix_pop_on;
    int burnin;
    int samples;
    int thinning;
    
    int rungs;
    std::vector<double> betaVec;
    std::vector<chain_noAdmixture> chainVec;
    
    chain_noAdmixture spareChain;
    std::vector<double> acceptanceRate;
    
    // likelihoods
    double logLikeGroup_sum;
    double logLikeGroup_sumSquared;
    double logLikeJoint_sum;
    double logLikeJoint_sumSquared;
    double harmonic;
    
    // Qmatrices
    std::vector< std::vector<double> > Qmatrix_ind;
    std::vector< std::vector<double> > Qmatrix_pop;
    
    std::vector< std::vector<double> > Qmatrix_spare;
    std::vector< std::vector<double> > Qmatrix_spare2;
    
    
    // PUBLIC FUNCTIONS
    
    // constructor
    MCMC_main_noAdmixture(globals &globals, int _Kindex, int _rungs);
    void reset (bool reset_Qmatrix_running);
    
    // perform MCMC
    void perform_MCMC(globals &globals, bool drawAlleleFreqs, bool fixLabels, bool outputLikelihood, bool outputPosteriorGrouping, int mainRep);
    

};

#endif
