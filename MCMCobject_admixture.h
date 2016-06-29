//
//  MavericK
//  MCMCobject_admixture.h
//
//  Created: Bob on 06/11/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Defines a class that can be used to carry out MCMC under the admixture model. The admixture parameter alpha can be defined as fixed or free to vary.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__MCMCobject_admixture__
#define __Maverick1_0__MCMCobject_admixture__

#include <iostream>
#include "globals.h"
#include "probability.h"
#include "misc.h"
#include "Hungarian.h"

//------------------------------------------------
// class containing all elements required for MCMC under admixture model
class MCMCobject_admixture {
    
public:
    
    // PUBLIC OBJECTS
    
    // basic quantities (copied over from globals object)
    std::vector< std::vector< std::vector<int> > > data;
    
    int Kindex;
    int K;
    int n;
    int loci;
    std::vector<int> J;
    std::vector<int> ploidy_vec;
    std::vector<std::string> uniquePops;
    int geneCopies;
    
    double lambda;
    bool fixAlpha_on;
    double alpha;
    double alphaPropSD;
    double beta;
    
    bool outputQmatrix_pop_on;
    
    int burnin;
    int samples;
    int thinning;
    
    std::vector< std::vector<double> > log_lookup;
    
    std::vector<int> linearGroup;
    std::vector< std::vector< std::vector<int> > > group;
    int groupIndex;
    std::vector< std::vector< std::vector<int> > > alleleCounts;
    std::vector< std::vector<int> > alleleCountsTotals;
    std::vector< std::vector< std::vector<double> > > alleleFreqs;
    std::vector< std::vector< std::vector<int> > > old_alleleCounts;
    std::vector< std::vector<int> > old_alleleCountsTotals;
    
    std::vector< std::vector<int> > admixCounts;
    std::vector<int> admixCountsTotals;
    std::vector< std::vector<double> > admixFreqs;
    std::vector< std::vector<int> > old_admixCounts;
    
    // likelihoods
    double logLikeGroup;
    double logLikeGroup_sum;
    std::vector<double> logLikeGroup_store;
    double logLikeGroup_sumSquared;
    double logLikeJoint;
    double logLikeJoint_sum;
    double logLikeJoint_sumSquared;
    double harmonic;
    
    std::vector<double> logProbVec;
    double logProbVecSum;
    double logProbVecMax;
    std::vector<double> probVec;
    double probVecSum;
    
    // Qmatrices. logQmatrix_ind_old, logQmatrix_ind_new and logQmatrix_running are used throughout MCMC (including burn-in phase) when solving label switching problem. Other Qmatrix objects are final outputs, and are only produced after burn-in phase.
    std::vector< std::vector<double> > logQmatrix_gene_old;
    std::vector< std::vector<double> > logQmatrix_gene_new;
    std::vector< std::vector<double> > Qmatrix_gene_new;
    std::vector< std::vector<double> > logQmatrix_gene_running;
    
    std::vector< std::vector<double> > logQmatrix_gene;
    std::vector< std::vector<double> > Qmatrix_gene;
    std::vector< std::vector<double> > Qmatrix_ind;
    std::vector< std::vector<double> > Qmatrix_pop;
    
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
    MCMCobject_admixture(globals &globals, int _Kindex, int _burnin, int _samples, int _thinning, double _beta);
    
    // perform MCMC
    void reset(bool reset_Qmatrix_running);
    void perform_MCMC(globals &globals, bool drawAlleleFreqs, bool storeLoglike, bool fixLabels, bool outputLikelihood, bool outputPosteriorGrouping, int mainRep);
    
    // update objects
    void group_update();
    void group_update_indLevel();
    void drawFreqs();
    void alpha_update();
    
    // label switching
    void chooseBestLabelPermutation(globals &globals, int rep);
    void produceQmatrix();
    void updateQmatrix(int rep);
    void storeQmatrix();
    
    // likelihoods
    void d_logLikeConditional(int i, int k, int targetGroup);
    void d_logLikeGroup();
    void d_logLikeJoint();
    
};

#endif
