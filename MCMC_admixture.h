//
//  MavericK
//  MCMC_admixture.h
//
//  Created: Bob on 06/11/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Defines a class that can be used to carry out MCMC under the admixture model. The admixture parameter alpha can be defined as fixed or free to vary.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__MCMC_admixture__
#define __Maverick1_0__MCMC_admixture__

#include <iostream>
#include "globals.h"
#include "probability.h"
#include "misc.h"
#include "Hungarian.h"
#include "particle_admixture.h"

//------------------------------------------------
// class containing all elements required for MCMC under admixture model
class MCMC_admixture {
    
public:
    
    // PUBLIC OBJECTS
    
    // basic quantities
    int Kindex;
    int K;
    int n;
    int loci;
    std::vector<int> ploidy_vec;
    bool outputQmatrix_pop_on;
    bool outputLikelihood_on;
    int nPops;
    int burnin;
    int samples;
    double GTI_pow;
    
    int rungs;
    std::vector<int> rung_order;
    std::vector<double> beta_raised_vec;
    std::vector<particle_admixture> particle_vec;
    std::vector<double> coupling_accept;
    
    // likelihoods for each rung
    std::vector<double> logLikeGroup_sum;
    std::vector<double> logLikeGroup_sumSquared;
    std::vector< std::vector<double> > logLikeGroup_store;
    double logLikeJoint_sum;
    double logLikeJoint_sumSquared;
    
    // Qmatrices
    std::vector< std::vector< std::vector< std::vector<double> > > > logQmatrix_gene_running;
    std::vector< std::vector< std::vector< std::vector<double> > > > logQmatrix_gene_update;
    std::vector< std::vector< std::vector< std::vector<double> > > > Qmatrix_gene_update;
    
    std::vector< std::vector< std::vector< std::vector<double> > > > Qmatrix_gene;
    std::vector< std::vector<double> > Qmatrix_ind;
    std::vector< std::vector<double> > Qmatrix_pop;
    
    // ordering of labels
    std::vector<int> labelOrder;
    std::vector<int> labelOrder_new;
    
    // autocorrelation
    std::vector<double> autoCorr;
    std::vector<double> ESS;
    
    // TI point estimates for each rung
    std::vector<double> TIpoint_mean;
    std::vector<double> TIpoint_var;
    std::vector<double> TIpoint_SE;
    
    // overall TI estimate and SE
    std::vector<double> logEvidence_TI_store;
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
    MCMC_admixture(globals &globals, int _Kindex);
    
    // perform MCMC
    void perform_MCMC(globals &globals);
    void MetropolisCoupling();
    void updateQmatrix(particle_admixture &particle, bool outOfBurnin);
    
};

#endif
