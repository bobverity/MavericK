//
//  MavericK
//  particle_admixture.h
//
//  Created: Bob on 22/06/2016
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Defines a class that can be used to carry out MCMC under the without-admixture model for a single value of beta. This class is a "particle" in the sense that it moves around updating it's various parameter values without retaining any information. All saving of parameter values, likelihoods etc. must occur at a higher level.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__particle_admixture__
#define __Maverick1_0__particle_admixture__

#include <iostream>
#include "globals.h"
#include "probability.h"
#include "misc.h"
#include "Hungarian.h"

//------------------------------------------------
// class containing all elements required for MCMC under no-admixture model
class particle_admixture {
    
public:
    
    // PUBLIC OBJECTS
    
    // basic quantities
    std::vector< std::vector< std::vector<int> > > data;
    int thisData;
    int K;
    int n;
    int loci;
    std::vector<int> J;
    std::vector<int> ploidy_vec;
    double lambda;
    double alpha;
    double beta;
    
    // grouping
    std::vector< std::vector< std::vector<int> > > group;
    int thisGroup;
    
    // allele and admix counts
    std::vector< std::vector< std::vector<int> > > alleleCounts;
    std::vector< std::vector<int> > alleleCountsTotals;
    std::vector< std::vector< std::vector<double> > > alleleFreqs;
    
    std::vector< std::vector< std::vector<int> > > old_alleleCounts;
    std::vector< std::vector<int> > old_alleleCountsTotals;
    
    std::vector< std::vector<int> > admixCounts;
    std::vector<int> admixCountsTotals;
    std::vector< std::vector<double> > admixFreqs;
    std::vector< std::vector<int> > old_admixCounts;
    
    // assignment probabilities
    std::vector<double> logProbVec;
    double logProbVecSum;
    double logProbVecMax;
    std::vector<double> probVec;
    double probVecSum;
    
    // likelihoods
    double logLikeGroup;
    double logLikeJoint;
    
    // Qmatrices. logQmatrix_gene_old, logQmatrix_gene_new and logQmatrix_gene_running are used throughout MCMC (including burn-in phase) when solving label switching problem. Other Qmatrix objects are final outputs, and are only produced after burn-in phase.
    std::vector< std::vector< std::vector< std::vector<double> > > > logQmatrix_gene_old;
    std::vector< std::vector< std::vector< std::vector<double> > > > logQmatrix_gene_new;
    std::vector< std::vector< std::vector< std::vector<double> > > > Qmatrix_gene_new;
    std::vector< std::vector< std::vector< std::vector<double> > > > logQmatrix_gene_running;
    
    std::vector< std::vector< std::vector< std::vector<double> > > > logQmatrix_gene;
    std::vector< std::vector< std::vector< std::vector<double> > > > Qmatrix_gene;
    
    // objects for Hungarian algorithm
    std::vector< std::vector<double> > costMat;
    std::vector<int> bestPerm;
    std::vector<int> bestPermOrder;
    
    std::vector<int>edgesLeft;
    std::vector<int>edgesRight;
    std::vector<int>blockedLeft;
    std::vector<int>blockedRight;
    
    // PUBLIC FUNCTIONS
    
    // constructors
    particle_admixture();
    particle_admixture(globals &globals, int _K, double _alpha, double _beta);
    
    // reset
    void reset(bool reset_Qmatrix_running);
    
    // update objects
    void group_update();
    void drawFreqs();
    
    // label switching
    void chooseBestLabelPermutation(globals &globals);
    void updateQmatrix();
    void storeQmatrix();
    
    // likelihoods
    void d_logLikeConditional(int i, int k);
    void d_logLikeGroup();
    void d_logLikeJoint();

};

#endif
