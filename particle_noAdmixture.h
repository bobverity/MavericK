//
//  MavericK
//  particle_noAdmixture.h
//
//  Created: Bob on 22/06/2016
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Defines a class that can be used to carry out MCMC under the without-admixture model for a single value of beta. This class is a "particle" in the sense that it moves around updating it's various parameter values without retaining any information. All saving of parameter values, likelihoods etc. must occur at a higher level.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__particle_noAdmixture__
#define __Maverick1_0__particle_noAdmixture__

#include <iostream>
#include "globals.h"
#include "probability.h"
#include "misc.h"
#include "Hungarian.h"

//------------------------------------------------
// class containing all elements required for MCMC under no-admixture model
class particle_noAdmixture {
    
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
    double beta_raised;
    
    // grouping
    std::vector<int> group;
    std::vector<int> group_propose;
    std::vector<int> group_order;
    
    // allele counts and frequencies
    std::vector< std::vector< std::vector<int> > > alleleCounts;
    std::vector< std::vector<int> > alleleCountsTotals;
    std::vector< std::vector< std::vector<double> > > alleleFreqs;
    
    // assignment probabilities
    std::vector<double> logProbVec;
    double logProbVecMax;
    std::vector<double> probVec;
    double probVecSum;
    
    // likelihoods
    double logLikeGroup;
    double logLikeJoint;
    
    // Qmatrices
    std::vector< std::vector<double> > Qmatrix_ind;
    std::vector< std::vector<double> > logQmatrix_ind;
    
    // PUBLIC FUNCTIONS
    
    // constructors
    particle_noAdmixture();
    particle_noAdmixture(globals &globals, int _K, double _beta_raised);
    
    ~particle_noAdmixture();
    
    // reset
    void reset();
    
    // update objects
    void group_update();
    void group_probs(int ind);
    void group_update_Klevel();
    void drawFreqs();
    
    // likelihoods
    void d_logLikeConditional(int i, int k);
    void d_logLikeGroup();
    void d_logLikeJoint();
    
};

#endif
