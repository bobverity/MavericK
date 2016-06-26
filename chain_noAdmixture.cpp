//
//  MavericK
//  chain_noAdmixture.cpp
//
//  Created: Bob on 22/06/2016
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "chain_noAdmixture.h"

using namespace std;

extern vector< vector<double> > log_lookup;
extern vector<double> log_lookup_0;

//------------------------------------------------
// chain_noAdmixture::
// default constructor for class
chain_noAdmixture::chain_noAdmixture(){};

//------------------------------------------------
// chain_noAdmixture::
// informed constructor for class
chain_noAdmixture::chain_noAdmixture(globals &globals, int _K, double _beta) {
    
    // copy some values over from globals object and arguments
    data = globals.data;
    K = _K;
    n = globals.n;
    loci = globals.loci;
    J = globals.J;
    ploidy_vec = globals.ploidy_vec;
    lambda = globals.lambda;
    beta = _beta;
    
    // initialise grouping vector
    group = vector<int>(n,1);
    
    // initialise allele counts and frequencies
    alleleCounts = vector< vector< vector<int> > >(K);
    alleleCountsTotals = vector< vector<int> >(K);
    alleleFreqs = vector< vector< vector<double> > >(K);
    for (int k=0; k<K; k++) {
        alleleCounts[k] = vector< vector<int> >(loci);
        alleleCountsTotals[k] = vector<int>(loci);
        alleleFreqs[k] = vector< vector<double> >(loci);
        for (int l=0; l<loci; l++) {
            alleleCounts[k][l] = vector<int>(J[l]);
            alleleFreqs[k][l] = vector<double>(J[l]);
        }
    }
    
    // initialise objects for calculating assignment probabilities
    logProbVec = vector<double>(K);
    logProbVecSum = 0;  // (used in Qmatrix calculation)
    logProbVecMax = 0;
    probVec = vector<double>(K);
    probVecSum = 0;
    
    // initialise Qmatrices
    logQmatrix_ind_old = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_ind_new = vector< vector<double> >(n, vector<double>(K));
    logQmatrix_ind_new = vector< vector<double> >(n, vector<double>(K));
    logQmatrix_ind_running = vector< vector<double> >(n, vector<double>(K));
    
    logQmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    
    // initialise objects for Hungarian algorithm
    costMat = vector< vector<double> >(K, vector<double>(K));
    bestPermOrder = vector<int>(K);
    
    edgesLeft = vector<int>(K);
    edgesRight = vector<int>(K);
    blockedLeft = vector<int>(K);
    blockedRight = vector<int>(K);
    
}

//------------------------------------------------
// chain_noAdmixture::
// reset objects used in MCMC
void chain_noAdmixture::reset(bool reset_Qmatrix_running) {
    
    // initialise group with random allocation
    vector<double> equalK(K,1/double(K));
    for (int i=0; i<n; i++) {
        group[i] = sample1(equalK,1.0);
    }
    
    // zero allele counts
    for (int k=0; k<K; k++) {
        alleleCountsTotals[k] = vector<int>(loci);
        for (int l=0; l<loci; l++) {
            alleleCounts[k][l] = vector<int>(J[l]);
        }
    }
    
    // populate allele counts based on grouping
    for (int ind=0; ind<n; ind++) {
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                if (data[ind][l][p]!=0) {
                    alleleCounts[group[ind]-1][l][data[ind][l][p]-1]++;
                    alleleCountsTotals[group[ind]-1][l]++;
                }
            }
        }
    }
    
    // reset likelihoods
    logLikeGroup = 0;
    logLikeJoint = 0;
    
    // reset Qmatrices
    logQmatrix_ind_old = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_ind_new = vector< vector<double> >(n, vector<double>(K));
    logQmatrix_ind_new = vector< vector<double> >(n, vector<double>(K));
    if (reset_Qmatrix_running) {
        logQmatrix_ind_running = vector< vector<double> >(n, vector<double>(K,-log(double(K))));
    }
    
    logQmatrix_ind = vector< vector<double> >(n, vector<double>(K, log(double(0))));
    Qmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    
}

//------------------------------------------------
// chain_noAdmixture::
// resample group allocation of all individuals by drawing from conditional posterior
void chain_noAdmixture::group_update() {
    
    // update group allocation for all individuals
    for (int ind=0; ind<n; ind++) {
        
        // subtract individual ind from allele counts
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                if (data[ind][l][p]!=0) {   // if not missing data
                    alleleCounts[group[ind]-1][l][data[ind][l][p]-1]--;
                    alleleCountsTotals[group[ind]-1][l]--;
                }
            }
        }
        
        // calculate probability of individual ind from all demes
        if (beta==0) {    // special case if beta==0 (draw from prior)
            logProbVec = vector<double>(K,-log(double(K)));
            probVec = vector<double>(K,1/double(K));
            probVecSum = 1.0;
        } else {
            for (int k=0; k<K; k++) {
                d_logLikeConditional(ind, k);   // update logProbVec[k]
            }
            logProbVecMax = *max_element(begin(logProbVec),end(logProbVec));
            probVecSum = 0;
            for (int k=0; k<K; k++) {
                probVec[k] = exp(beta*logProbVec[k]-beta*logProbVecMax);
                probVecSum += probVec[k];
            }
            // use these values to define Qmatrix_ind_new for this iteration
            logProbVecSum = log(probVecSum);
            for (int k=0; k<K; k++) {
                logQmatrix_ind_new[ind][k] = logProbVec[k]- (logProbVecSum + logProbVecMax);
                Qmatrix_ind_new[ind][k] = probVec[k]/probVecSum;
            }
        }
        
        // resample grouping
        group[ind] = sample1(probVec, probVecSum);
        
        // add individual ind to allele counts
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                if (data[ind][l][p]!=0) {   // if not missing data
                    alleleCounts[group[ind]-1][l][data[ind][l][p]-1]++;
                    alleleCountsTotals[group[ind]-1][l]++;
                }
            }
        }
        
    }
    
}

//------------------------------------------------
// chain_noAdmixture::
// calculate logQmatrix_ind_new for this iteration (commented out because using short-cut method of fixing label switching)
void chain_noAdmixture::produceQmatrix() {
    /*
    // populate Qmatrix_ind_new
    for (int i=0; i<n; i++) {
        logProbVecSum = log(double(0));
        for (int k=0; k<K; k++) {
            d_logLikeConditional(i, k);   // update logProbVec[k]
            logProbVecSum = logSum(logProbVecSum, logProbVec[k]);
        }
        for (int k=0; k<K; k++) {
            logQmatrix_ind_new[i][k] = logProbVec[k]-logProbVecSum;
        }
    }
    */
}

//------------------------------------------------
// chain_noAdmixture::
// choose best permutation of labels using method of Stephens (2000)
void chain_noAdmixture::chooseBestLabelPermutation(globals &globals) {
    
    // calculate cost matrix from old and new Qmatrices
    for (int k1=0; k1<K; k1++) {
        for (int k2=0; k2<K; k2++) {
            costMat[k1][k2] = 0;
            for (int i=0; i<n; i++) {
                costMat[k1][k2] += Qmatrix_ind_new[i][k1]*(logQmatrix_ind_new[i][k1]-logQmatrix_ind_running[i][k2]);
            }
        }
    }
    
    // find best permutation of current labels
    bestPerm = hungarian(costMat, edgesLeft, edgesRight, blockedLeft, blockedRight, globals.outputLog_on, globals.outputLog_fileStream);

    // define bestPermOrder. If the numbers 1:m_K are placed in best-perm-order then we arrive back at bestPerm. In R terms we would say bestPermOrder=order(bestPerm).
    bool performSwap = false;
    for (int k=0; k<K; k++) {
        bestPermOrder[bestPerm[k]] = k;
        if (bestPerm[k]!=k)
            performSwap = true;
    }
    
    // swap labels if necessary
    if (performSwap) {
        
        // update grouping to reflect swapped labels
        for (int i=0; i<n; i++) {
            group[i] = bestPerm[group[i]-1]+1;
        }
        
        // update allele counts to reflect swapped labels
        old_alleleCounts = alleleCounts;
        old_alleleCountsTotals = alleleCountsTotals;
        for (int k=0; k<K; k++) {
            alleleCounts[k] = old_alleleCounts[bestPermOrder[k]];
            alleleCountsTotals[k] = old_alleleCountsTotals[bestPermOrder[k]];
        }
        
        // update logQmatrix_ind_new to reflect swapped labels
        logQmatrix_ind_old = logQmatrix_ind_new;
        for (int i=0; i<n; i++) {
            for (int k=0; k<K; k++) {
                logQmatrix_ind_new[i][k] = logQmatrix_ind_old[i][bestPermOrder[k]];
            }
        }
        
    }
    
}

//------------------------------------------------
// chain_noAdmixture::
// add logQmatrix_ind_new to logQmatrix_ind_running
void chain_noAdmixture::updateQmatrix() {
    
    for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
            logQmatrix_ind_running[i][k] = logSum(logQmatrix_ind_running[i][k], logQmatrix_ind_new[i][k]);
        }
    }
}

//------------------------------------------------
// chain_noAdmixture::
// store Qmatrix values
void chain_noAdmixture::storeQmatrix() {
    
    // store individual-level Qmatrix
    for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
            logQmatrix_ind[i][k] = logSum(logQmatrix_ind[i][k], logQmatrix_ind_new[i][k]);
        }
    }
    
}

//------------------------------------------------
// chain_noAdmixture::
// probability of data given grouping only, integrated over unknown allele frequencies
void chain_noAdmixture::d_logLikeGroup() {
    
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
// chain_noAdmixture::
// draw allele frequencies given allele counts and lambda prior
void chain_noAdmixture::drawFreqs() {
    
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
    
}

//------------------------------------------------
// chain_noAdmixture::
// probability of data given grouping and allele frequencies
void chain_noAdmixture::d_logLikeJoint() {
    
    // calculate likelihood
    logLikeJoint = 0;
    double running = 1.0;
    for (int i=0; i<n; i++) {
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[i]; p++) {
                if (data[i][l][p]!=0) {
                    running *= alleleFreqs[group[i]-1][l][data[i][l][p]-1];
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

//------------------------------------------------
// chain_noAdmixture::
// conditional probability of ith individual from kth deme (output in log space)
void chain_noAdmixture::d_logLikeConditional(int i, int k) {
    
    // calculate conditional probability of data
    logProbVec[k] = 0;
    int d, a, a_t;  // for making temporary copies of data, alleleCounts, and alleleCountsTotals respectively
    double log_part2;
    for (unsigned int l=0; l<loci; l++) {
        a_t = alleleCountsTotals[k][l];
        log_part2 = log_lookup[a_t][J[l]-1];
        for (unsigned int p=0; p<ploidy_vec[i]; p++) {
            d = data[i][l][p];
            a = alleleCounts[k][l][d-1];
            if (d!=0) {
                if ((a<int(1e4)) && (a_t<int(1e4))) {
                    logProbVec[k] += log_lookup_0[a]-log_part2;
                } else {
                    logProbVec[k] += log((a + lambda)/double(a_t + J[l]*lambda));
                }
                alleleCounts[k][l][d-1] ++;
                a_t ++;
            }
        }
        for (unsigned int p=0; p<ploidy_vec[i]; p++) {
            d = data[i][l][p];
            if (d!=0) {
                alleleCounts[k][l][d-1] --;
            }
        }
    }
}



