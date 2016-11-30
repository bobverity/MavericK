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
    bool outputQmatrix_pop_on;
    int nPops;
    int burnin;
    int samples;
    
    int rungs;
    std::vector<int> rungOrder;
    int rung1;
    int rung2;
    std::vector<double> betaVec;
    std::vector<particle_noAdmixture> particleVec;
    std::vector<double> acceptanceRate;
    
    // likelihoods for each rung
    std::vector<double> logLikeGroup_sum;
    std::vector<double> logLikeGroup_sumSquared;
    std::vector< std::vector<double> > logLikeGroup_store;
    double logLikeJoint_sum;
    double logLikeJoint_sumSquared;
    
    // Qmatrices
    std::vector< std::vector<double> > logQmatrix_ind_running;
    std::vector< std::vector<double> > logQmatrix_ind_update;
    std::vector< std::vector<double> > Qmatrix_ind_update;
    
    std::vector< std::vector<double> > Qmatrix_ind;
    std::vector< std::vector<double> > Qmatrix_pop;
    
    // ordering of labels
    std::vector<int> labelOrder;
    std::vector<int> labelOrder_new;
    
    // autocorrelation
    std::vector<double> autoCorr;
    std::vector<double> ESS;
    
    // harmonic mean estimator
    double harmonic;
    
    // Structure estimator
    double structure_loglike_mean;
    double structure_loglike_var;
    double logEvidence_structure;
    
    // TI point estimates for each rung
    std::vector<double> TIpoint_mean;
    std::vector<double> TIpoint_var;
    std::vector<double> TIpoint_SE;
    
    // overall TI point estimate and SE
    double logEvidence_TI;
    double logEvidence_TI_var;
    double logEvidence_TI_SE;
    
    // objects for Hungarian algorithm
    std::vector< std::vector<double> > costMat;
    std::vector<int> bestPerm;
    std::vector<int> bestPermOrder;
    
    std::vector<int>edgesLeft;
    std::vector<int>edgesRight;
    std::vector<int>blockedLeft;
    std::vector<int>blockedRight;
    
    
    // PUBLIC FUNCTIONS
    
    // constructor
    MCMC_noAdmixture(globals &globals, int _Kindex, int _rungs);
    
    // perform MCMC
    void perform_MCMC(globals &globals);
    void MetropolisCoupling();
    void updateQmatrix(particle_noAdmixture &particle, bool outOfBurnin);

};

#endif
