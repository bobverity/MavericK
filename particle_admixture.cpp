//
//  MavericK
//  particle_admixture.cpp
//
//  Created: Bob on 22/06/2016
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "particle_admixture.h"

using namespace std;

extern vector< vector<double> > log_lookup;
extern vector<double> log_lookup_0;

extern int test1;
extern int test2;

//------------------------------------------------
// particle_admixture::
// default constructor for class
particle_admixture::particle_admixture(){};

//------------------------------------------------
// particle_admixture::
// informed constructor for class
particle_admixture::particle_admixture(globals &globals, int _K, double _alpha, double _alphaPropSD, double _beta) {
    
    // copy some values over from globals object and arguments
    data = globals.data;
    K = _K;
    n = globals.n;
    loci = globals.loci;
    J = globals.J;
    ploidy_vec = globals.ploidy_vec;
    lambda = globals.lambda;
    beta = _beta;
    alpha = _alpha;
    alphaPropSD = _alphaPropSD;
    
    // initialise grouping array
    group = vector< vector< vector<int> > >(n, vector< vector<int> >(loci));
    group_order = vector<int>(n);
    for (int ind=0; ind<n; ind++) {
        group_order[ind] = ind+1;
        for (int l=0; l<loci; l++) {
            group[ind][l] = vector<int>(ploidy_vec[ind]);
        }
    }
    group_propose = group;
    
    // initialise allele counts and frequencies
    alleleCounts = vector< vector< vector<int> > >(K, vector< vector<int> >(loci));
    alleleCountsTotals = vector< vector<int> >(K, vector<int>(loci));
    alleleFreqs = vector< vector< vector<double> > >(K, vector< vector<double> >(loci));
    for (int k=0; k<K; k++) {
        for (int l=0; l<loci; l++) {
            alleleCounts[k][l] = vector<int>(J[l]);
            alleleFreqs[k][l] = vector<double>(J[l]);
        }
    }
    
    // initialise admix counts and frequencies
    admixCounts = vector< vector<int> >(n,vector<int>(K));
    admixCountsTotals = vector<int>(n);
    admixFreqs = vector< vector<double> >(n,vector<double>(K));
    
    // initialise objects for calculating assignment probabilities
    logProbVec = vector<double>(K);
    logProbVecMax = 0;
    probVec = vector<double>(K);
    probVecSum = 0;
    
    // initialise likelihoods
    logLikeGroup = 0;
    logLikeJoint = 0;
    
    // initialise Qmatrices
    logQmatrix_gene = vector< vector< vector< vector<double> > > >(n, vector< vector< vector<double> > >(loci));
    Qmatrix_gene = vector< vector< vector< vector<double> > > >(n, vector< vector< vector<double> > >(loci));
    for (int ind=0; ind<n; ind++) {
        logQmatrix_gene[ind] = vector< vector< vector<double> > >(loci, vector< vector<double> >(ploidy_vec[ind], vector<double>(K)));
        Qmatrix_gene[ind] = vector< vector< vector<double> > >(loci, vector< vector<double> >(ploidy_vec[ind], vector<double>(K)));
    }
    
}

//------------------------------------------------
// particle_admixture::
// reset objects used in MCMC
void particle_admixture::reset(bool reset_Qmatrix_running) {
    
    // initialise group with random allocation
    for (int ind=0; ind<n; ind++) {
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                group[ind][l][p] = sample2(1,K);
            }
        }
    }
    
    // zero allele counts and admix counts
    for (int k=0; k<K; k++) {
        alleleCountsTotals[k] = vector<int>(loci);
        for (int l=0; l<loci; l++) {
            for (int j=0; j<J[l]; j++) {
                alleleCounts[k][l][j] = 0;
            }
        }
    }
    for (int ind=0; ind<n; ind++) {
        admixCountsTotals[ind]= 0;
        for (int k=0; k<K; k++) {
            admixCounts[ind][k] = 0;
        }
    }
    
    // populate allele counts and admix counts based on grouping
    for (int ind=0; ind<n; ind++) {
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                thisData = data[ind][l][p];
                if (thisData!=0) {
                    thisGroup = group[ind][l][p];
                    alleleCounts[thisGroup-1][l][thisData-1]++;
                    alleleCountsTotals[thisGroup-1][l]++;
                    
                    admixCounts[ind][thisGroup-1]++;
                    admixCountsTotals[ind]++;
                }
            }
        }
    }
    
    // reset likelihoods
    logLikeGroup = 0;
    logLikeJoint = 0;
    
    // reset Qmatrices
    for (int ind=0; ind<n; ind++) {
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                for (int k=0; k<K; k++) {
                    logQmatrix_gene[ind][l][p][k] = 0;
                    Qmatrix_gene[ind][l][p][k] = 0;
                }
            }
        }
    }
    
}

//------------------------------------------------
// particle_admixture::
// resample group allocation of all individuals by drawing from conditional posterior
void particle_admixture::group_update() {
    
    // update group allocation for all gene copies
    for (int ind=0; ind<n; ind++) {
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                thisGroup = group[ind][l][p];
                thisData = data[ind][l][p];
                
                // subtract this gene copy from allele counts and admix counts
                if (thisData!=0) {   // if not missing data
                    alleleCounts[thisGroup-1][l][thisData-1]--;
                    alleleCountsTotals[thisGroup-1][l]--;
                    
                    admixCounts[ind][thisGroup-1]--;
                    admixCountsTotals[ind]--;
                }
                
                // resample group allocation of individual ind
                if (beta==0) {    // special case if beta==0 (draw from prior)
                    thisGroup = sample2(1,K);
                } else {
                    group_probs(ind,l,p);
                    thisGroup = sample1(probVec, probVecSum);
                }
                group[ind][l][p] = thisGroup;
                
                // add this gene copy to allele counts and admix counts
                if (thisData!=0) {   // if not missing data
                    alleleCounts[thisGroup-1][l][thisData-1]++;
                    alleleCountsTotals[thisGroup-1][l]++;
                    
                    admixCounts[ind][thisGroup-1]++;
                    admixCountsTotals[ind]++;
                }
                
            } // p
        } // l
    } // ind
    
}

//------------------------------------------------
// particle_admixture::
// calculate conditional probability of this gene copy coming from each deme
void particle_admixture::group_probs(int ind, int l, int p) {
    
    thisData = data[ind][l][p];
    probVecSum = 0;
    for (int k=0; k<K; k++) {
        if (thisData==0) {
            logProbVec[k] = 0;
        } else {
            thisAlleleCounts = alleleCounts[k][l][thisData-1];
            thisAlleleCountsTotals = alleleCountsTotals[k][l];
            if ((thisAlleleCounts<int(1e4)) && (thisAlleleCountsTotals<int(1e4))) {
                logProbVec[k] = beta*log_lookup_0[thisAlleleCounts] - log_lookup[thisAlleleCountsTotals][J[l]-1];
            } else {
                logProbVec[k] = beta*log((thisAlleleCounts + lambda)/double(thisAlleleCountsTotals + J[l]*lambda));
            }
        }
        probVec[k] = double(admixCounts[ind][k]+alpha)*exp(logProbVec[k]);
        probVecSum += probVec[k];
    }
    
}

//------------------------------------------------
// particle_admixture::
// Metropolis-Hastings step to update all gene copies within an individual simultaneously. This helps with mixing when alpha is small, as otherwise it can be very difficult for an individual allocated to the wrong group to move freely.
void particle_admixture::group_update_indLevel() {
    
    double logLike_old;
    double logLike_new;
    double propose_logProb_old;
    double propose_logProb_new;
    double MH_diff;
    double rand1;
    
    // loop over all individuals
    for (int ind=0; ind<n; ind++) {
        
        // subtract all gene copies in this individual and calculate likelihood of current grouping at the same time
        logLike_old = 0;
        propose_logProb_old = 0;
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                thisData = data[ind][l][p];
                thisGroup = group[ind][l][p];
                
                // subtract this gene copy from allele counts and admix counts
                if (thisData!=0) {   // if not missing data
                    alleleCounts[thisGroup-1][l][thisData-1]--;
                    alleleCountsTotals[thisGroup-1][l]--;
                    
                    admixCounts[ind][thisGroup-1]--;
                    admixCountsTotals[ind]--;
                }
                
                // calculate probability of this gene copy from all demes
                group_probs(ind,l,p);
                
                // calculate probability of chosen grouping
                propose_logProb_old += log(probVec[thisGroup-1]) - log(probVecSum);
                logLike_old += log(probVec[thisGroup-1]);
                
            }
        }
        
        // propose a new grouping and calculate likelihood
        logLike_new = 0;
        propose_logProb_new = 0;
        for (unsigned int l=0; l<loci; l++) {
            for (unsigned int p=0; p<ploidy_vec[ind]; p++) {
                thisData = data[ind][l][p];
                
                // calculate probability of this gene copy from all demes
                group_probs(ind,l,p);
                
                // resample grouping
                thisGroup = sample1(probVec, probVecSum);
                group_propose[ind][l][p] = thisGroup;
                
                // calculate probability of new grouping
                propose_logProb_new += log(probVec[thisGroup-1]) - log(probVecSum);
                logLike_new += log(probVec[thisGroup-1]);
                
                // add this gene copy to allele counts and admix counts
                if (thisData!=0) {   // if not missing data
                    alleleCounts[thisGroup-1][l][thisData-1]++;
                    alleleCountsTotals[thisGroup-1][l]++;
                    
                    admixCounts[ind][thisGroup-1]++;
                    admixCountsTotals[ind]++;
                }
            }
        }
        
        // Metropolis-Hastings step. If accept then stick with new grouping, otherwise revert back
        MH_diff = (logLike_new - propose_logProb_new) - (logLike_old - propose_logProb_old);
        rand1 = runif1(0,1);
        if (log(rand1)<MH_diff) {
            
            // accept move
            for (unsigned int l=0; l<loci; l++) {
                for (unsigned int p=0; p<ploidy_vec[ind]; p++) {
                    group[ind][l][p] = group_propose[ind][l][p];
                }
            }
            
        } else {
            
            // reject move
            for (unsigned int l=0; l<loci; l++) {
                for (unsigned int p=0; p<ploidy_vec[ind]; p++) {
                    thisData = data[ind][l][p];
                    if (thisData!=0) {   // if not missing data
                        
                        // subtract new group
                        thisGroup = group_propose[ind][l][p];
                        alleleCounts[thisGroup-1][l][thisData-1]--;
                        alleleCountsTotals[thisGroup-1][l]--;
                        
                        admixCounts[ind][thisGroup-1]--;
                        admixCountsTotals[ind]--;
                        
                        // reinstate old group
                        thisGroup = group[ind][l][p];
                        alleleCounts[thisGroup-1][l][thisData-1]++;
                        alleleCountsTotals[thisGroup-1][l]++;
                        
                        admixCounts[ind][thisGroup-1]++;
                        admixCountsTotals[ind]++;
                    }
                    
                }
            }
            
        }
        
    }   // end loop over individuals
    
}

//------------------------------------------------
// particle_admixture::
// resample group allocation of all individuals by drawing from conditional posterior
void particle_admixture::group_update_Klevel() {
    
    // if beta==0 then all moves accepted, so skip
    if (beta==0)
        return;
    
    // choose two distinct random groups
    int K1 = sample2(1,K);
    int K2 = sample2(1,K);
    if (K2==K1) {
        K2++;
        if (K2>K)
            K2 = 1;
    }
    
    // search through observations in random order
    shuffle1(group_order);
    
    // subtract all individuals in these two groups
    double logLike_old = 0;
    double propose_logProb_old = 0;
    for (int i=(n-1); i>=0; i--) {  // loop backwards through group_order
        int ind = group_order[i]-1;
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                thisGroup = group[ind][l][p];
                if (thisGroup==K1 || thisGroup==K2) {
                    
                    thisData = data[ind][l][p];
                    if (thisData==0) {
                        probVec[thisGroup-1] = 0.5;
                        probVecSum = 1;
                    } else {
                        alleleCounts[thisGroup-1][l][thisData-1]--;
                        alleleCountsTotals[thisGroup-1][l]--;
                        
                        admixCounts[ind][thisGroup-1]--;
                        admixCountsTotals[ind]--;
                        
                        thisAlleleCounts = alleleCounts[K1-1][l][thisData-1];
                        thisAlleleCountsTotals = alleleCountsTotals[K1-1][l];
                        if ((thisAlleleCounts<int(1e4)) && (thisAlleleCountsTotals<int(1e4))) {
                            logProbVec[K1-1] = beta*log_lookup_0[thisAlleleCounts] - log_lookup[thisAlleleCountsTotals][J[l]-1];
                        } else {
                            logProbVec[K1-1] = beta*log((thisAlleleCounts + lambda)/double(thisAlleleCountsTotals + J[l]*lambda));
                        }
                        probVec[K1-1] = double(admixCounts[ind][K1-1]+alpha)*exp(logProbVec[K1-1]);
                        
                        thisAlleleCounts = alleleCounts[K2-1][l][thisData-1];
                        thisAlleleCountsTotals = alleleCountsTotals[K2-1][l];
                        if ((thisAlleleCounts<int(1e4)) && (thisAlleleCountsTotals<int(1e4))) {
                            logProbVec[K2-1] = beta*log_lookup_0[thisAlleleCounts] - log_lookup[thisAlleleCountsTotals][J[l]-1];
                        } else {
                            logProbVec[K2-1] = beta*log((thisAlleleCounts + lambda)/double(thisAlleleCountsTotals + J[l]*lambda));
                        }
                        probVec[K2-1] = double(admixCounts[ind][K2-1]+alpha)*exp(logProbVec[K2-1]);
                        
                        probVecSum = probVec[K1-1] + probVec[K2-1];
                    }
                    
                    // calculate probability of chosen grouping
                    propose_logProb_old += log(probVec[thisGroup-1]/probVecSum);
                    logLike_old += log(probVec[thisGroup-1]);
                    
                }   // end if group==K1 or K2
            }   // end loop over ploidy
        }   // end loop over loci
    }   // end loop over individuals
    
    // propose a new grouping and calculate likelihood
    double logLike_new = 0;
    double propose_logProb_new = 0;
    for (int i=0; i<n; i++) {  // loop forwards through group_order
        int ind = group_order[i]-1;
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                thisGroup = group[ind][l][p];
                if (thisGroup==K1 || thisGroup==K2) {
                    
                    thisData = data[ind][l][p];
                    if (thisData==0) {
                        probVec[thisGroup-1] = 0.5;
                        probVecSum = 1;
                    } else {
                        
                        thisAlleleCounts = alleleCounts[K1-1][l][thisData-1];
                        thisAlleleCountsTotals = alleleCountsTotals[K1-1][l];
                        if ((thisAlleleCounts<int(1e4)) && (thisAlleleCountsTotals<int(1e4))) {
                            logProbVec[K1-1] = beta*log_lookup_0[thisAlleleCounts] - log_lookup[thisAlleleCountsTotals][J[l]-1];
                        } else {
                            logProbVec[K1-1] = beta*log((thisAlleleCounts + lambda)/double(thisAlleleCountsTotals + J[l]*lambda));
                        }
                        probVec[K1-1] = double(admixCounts[ind][K1-1]+alpha)*exp(logProbVec[K1-1]);
                        
                        thisAlleleCounts = alleleCounts[K2-1][l][thisData-1];
                        thisAlleleCountsTotals = alleleCountsTotals[K2-1][l];
                        if ((thisAlleleCounts<int(1e4)) && (thisAlleleCountsTotals<int(1e4))) {
                            logProbVec[K2-1] = beta*log_lookup_0[thisAlleleCounts] - log_lookup[thisAlleleCountsTotals][J[l]-1];
                        } else {
                            logProbVec[K2-1] = beta*log((thisAlleleCounts + lambda)/double(thisAlleleCountsTotals + J[l]*lambda));
                        }
                        probVec[K2-1] = double(admixCounts[ind][K2-1]+alpha)*exp(logProbVec[K2-1]);
                        
                        probVecSum = probVec[K1-1] + probVec[K2-1];
                        
                        // resample group allocation of individual ind
                        if (rbernoulli1(probVec[K1-1]/probVecSum)) {
                            group_propose[ind][l][p] = K1;
                        } else {
                            group_propose[ind][l][p] = K2;
                        }
                        
                        alleleCounts[group_propose[ind][l][p]-1][l][thisData-1]++;
                        alleleCountsTotals[group_propose[ind][l][p]-1][l]++;
                        
                        admixCounts[ind][group_propose[ind][l][p]-1]++;
                        admixCountsTotals[ind]++;
                    }
                    
                    // calculate probability of chosen grouping
                    propose_logProb_new += log(probVec[group_propose[ind][l][p]-1]/probVecSum);
                    logLike_new += log(probVec[group_propose[ind][l][p]-1]);
                    
                }   // end if group==K1 or K2
            }   // end loop over ploidy
        }   // end loop over loci
    }   // end loop over individuals
    
    // Metropolis-Hastings step. If accept then stick with new grouping, otherwise revert back
    double MH_diff = (logLike_new - propose_logProb_new) - (logLike_old - propose_logProb_old);
    double rand1 = runif1(0,1);
    
    test2++;
    
    if (log(rand1)<MH_diff) {
        
        test1++;
        
        // accept move
        for (int ind=0; ind<n; ind++) {
            for (int l=0; l<loci; l++) {
                for (int p=0; p<ploidy_vec[ind]; p++) {
                    thisGroup = group[ind][l][p];
                    if (thisGroup==K1 || thisGroup==K2) {
                        group[ind][l][p] = group_propose[ind][l][p];
                    }
                }
            }
        }
        
    } else {
        
        // reject move
        for (int ind=0; ind<n; ind++) {
            for (int l=0; l<loci; l++) {
                for (int p=0; p<ploidy_vec[ind]; p++) {
                    thisGroup = group[ind][l][p];
                    if (thisGroup==K1 || thisGroup==K2) {
                        thisData = data[ind][l][p];
                        if (thisData!=0) {   // if not missing data
                            
                            // subtract new group
                            alleleCounts[group_propose[ind][l][p]-1][l][thisData-1]--;
                            alleleCountsTotals[group_propose[ind][l][p]-1][l]--;
                            
                            admixCounts[ind][group_propose[ind][l][p]-1]--;
                            admixCountsTotals[ind]--;
                            
                            // reinstate old group
                            alleleCounts[thisGroup-1][l][thisData-1]++;
                            alleleCountsTotals[thisGroup-1][l]++;
                            
                            admixCounts[ind][thisGroup-1]++;
                            admixCountsTotals[ind]++;
                            
                        }
                    }
                }
            }
        }
        
    }
    
}

//------------------------------------------------
// particle_admixture::
// Gibbs sampler for alpha using data augmentation
void particle_admixture::alpha_update() {
    
    // prior parameters on alpha (shape and rate of gamma prior)
    double epsilon = 1.0;
    double lambda = 0.2;
    
    // draw latent variables
    int z = 0;
    double logPhi = 0;
    for (int i=0; i<n; i++) {
        logPhi += log(rbeta1(K*alpha+1, admixCountsTotals[i]));
        for (int k=0; k<K; k++) {
            z += rCRPgroups(admixCounts[i][k], alpha);
        }
    }
    
    // draw final gamma parameters
    double w1 = epsilon + z - n;
    double w2 = lambda - K*logPhi;
    double p;
    for (int i=0; i<n; i++) {
        p = admixCountsTotals[i]/(admixCountsTotals[i] + K*w1/w2);
        w1 += rbernoulli1(1-p);
    }
    
    // draw new alpha
    alpha = rgamma1(w1, w2);
}

//------------------------------------------------
// particle_admixture::
// probability of data given grouping only, integrated over unknown allele frequencies
void particle_admixture::d_logLikeGroup() {
    
    // Multinomial-Dirichlet likelihood
    logLikeGroup = 0;
    for (int k=0; k<K; k++) {
        for (int l=0; l<loci; l++) {
            for (int j=0; j<J[l]; j++) {
                logLikeGroup += lgamma(lambda + alleleCounts[k][l][j]) - lgamma(lambda);
            }
            logLikeGroup += lgamma(J[l]*lambda) - lgamma(J[l]*lambda + alleleCountsTotals[k][l]);
        }
    }
    
}

//------------------------------------------------
// particle_admixture::
// draw allele frequencies and admixture proportions, given allele and admixture counts and lambda and alpha values
void particle_admixture::drawFreqs() {
    
    // draw allele frequencies
    double randSum;
    for (int k=0; k<K; k++) {
        for (int l=0; l<loci; l++) {
            randSum = 0;
            for (int j=0; j<J[l]; j++) {
                alleleFreqs[k][l][j] = rgamma1(alleleCounts[k][l][j]+lambda, 1.0);
                randSum += alleleFreqs[k][l][j];
            }
            for (int j=0; j<J[l]; j++) {
                alleleFreqs[k][l][j] /= randSum;
            }
            
        }
    }
    
    // draw admixture proportions
    for (int ind=0; ind<n; ind++) {
        randSum = 0;
        for (int k=0; k<K; k++) {
            admixFreqs[ind][k] = rgamma1(admixCounts[ind][k]+alpha, 1.0);
            randSum += admixFreqs[ind][k];
        }
        for (int k=0; k<K; k++) {
            admixFreqs[ind][k] /= randSum;
        }
    }
    
}

//------------------------------------------------
// particle_admixture::
// probability of data given grouping and allele frequencies
void particle_admixture::d_logLikeJoint() {
    
    // calculate likelihood
    logLikeJoint = 0;
    double running = 1.0;
    for (int ind=0; ind<n; ind++) {
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                if (data[ind][l][p]!=0) {
                    running *= alleleFreqs[group[ind][l][p]-1][l][data[ind][l][p]-1];
                }
                if (running<UNDERFLO) {
                    logLikeJoint += log(running);
                    running = 1.0;
                }
            }
        }
    }
    logLikeJoint += log(running);
    
}




