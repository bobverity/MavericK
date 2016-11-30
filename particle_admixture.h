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
    int K;
    int n;
    int loci;
    std::vector<int> J;
    std::vector<int> ploidy_vec;
    double lambda;
    double beta;
    
    double alpha;
    double alphaPropSD;
    
    // grouping
    std::vector< std::vector< std::vector<int> > > group;
    std::vector< std::vector< std::vector<int> > > group_propose;
    std::vector<int> group_order;
    
    // allele counts and frequencies
    std::vector< std::vector< std::vector<int> > > alleleCounts;
    std::vector< std::vector<int> > alleleCountsTotals;
    std::vector< std::vector< std::vector<double> > > alleleFreqs;
    
    // admix counts and frequencies
    std::vector< std::vector<int> > admixCounts;
    std::vector<int> admixCountsTotals;
    std::vector< std::vector<double> > admixFreqs;
    
    // scalars for storing current status
    int thisData;
    int thisGroup;
    int thisAlleleCounts;
    int thisAlleleCountsTotals;
    
    // assignment probabilities
    std::vector<double> logProbVec;
    double logProbVecMax;
    std::vector<double> probVec;
    double probVecSum;
    
    // likelihoods
    double logLikeGroup;
    double logLikeJoint;
    
    // Qmatrices
    std::vector< std::vector< std::vector< std::vector<double> > > > logQmatrix_gene;
    std::vector< std::vector< std::vector< std::vector<double> > > > Qmatrix_gene;
    
    
    // PUBLIC FUNCTIONS
    
    // constructors
    particle_admixture();
    particle_admixture(globals &globals, int _K, double _alpha, double _alphaPropSD, double _beta);
    
    // reset
    void reset(bool reset_Qmatrix_running);
    
    // update objects
    void group_update();
    void group_probs(int ind, int l, int p);
    void group_update_indLevel();
    void group_update_Klevel();
    void alpha_update();
    void drawFreqs();
    
    // likelihoods
    void d_logLikeConditional(int i, int k);
    void d_logLikeGroup();
    void d_logLikeJoint();

};

#endif
