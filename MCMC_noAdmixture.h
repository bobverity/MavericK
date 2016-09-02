//
//  MavericK
//  MCMC_TI_noAdmixture.h
//
//  Created: Bob on 22/06/2016
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Defines a class that can be used to carry out MCMC under the without-admixture model.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__MCMC_noAdmixture__
#define __Maverick1_0__MCMC_noAdmixture__

#include <iostream>
#include "globals.h"
#include "probability.h"
#include "misc.h"
#include "Hungarian.h"
#include "particle_noAdmixture.h"

//------------------------------------------------
// class containing all elements required for MCMC under no-admixture model
class MCMC_noAdmixture {
    
public:
    
    // PUBLIC OBJECTS
    
    // basic quantities
    int Kindex;
    int K;
    int n;
    std::vector<std::string> uniquePops;
    bool outputQmatrix_pop_on;
    int burnin;
    int samples;
    
    int rungs;
    std::vector<double> betaVec;
    std::vector<particle_noAdmixture> particleVec;
    
    particle_noAdmixture spareParticle;
    std::vector<double> acceptanceRate;
    
    // likelihoods
    std::vector<double> logLikeGroup_sum;
    std::vector<double> logLikeGroup_sumSquared;
    std::vector< std::vector<double> > logLikeGroup_store;
    double logLikeJoint_sum;
    double logLikeJoint_sumSquared;
    double harmonic;
    
    // Qmatrices
    std::vector< std::vector<double> > Qmatrix_ind;
    std::vector< std::vector<double> > Qmatrix_pop;
    
    std::vector< std::vector<double> > Qmatrix_spare;
    std::vector< std::vector<double> > Qmatrix_spare2;
    
    // autocorrelation
    std::vector<double> autoCorr;
    std::vector<double> ESS;
    
    // TI point estimates for each rung
    std::vector<double> TIpoint_mean;
    std::vector<double> TIpoint_var;
    std::vector<double> TIpoint_SE;
    
    // overall TI point estimate and SE
    double logEvidence_TI;
    double logEvidence_TI_var;
    double logEvidence_TI_SE;
    
    // PUBLIC FUNCTIONS
    
    // constructor
    MCMC_noAdmixture(globals &globals, int _Kindex, int _rungs);
    
    // perform MCMC
    void perform_MCMC(globals &globals, bool outputLikelihood, bool outputPosteriorGrouping);
    

};

#endif
