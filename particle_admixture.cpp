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

//------------------------------------------------
// particle_admixture::
// default constructor for class
particle_admixture::particle_admixture(){};

//------------------------------------------------
// particle_admixture::
// informed constructor for class
particle_admixture::particle_admixture(globals &globals, int _K, double _alpha, double _beta) {
    
    // copy some values over from globals object and arguments
    data = globals.data;
    K = _K;
    n = globals.n;
    loci = globals.loci;
    J = globals.J;
    ploidy_vec = globals.ploidy_vec;
    lambda = globals.lambda;
    alpha = _alpha;
    beta = _beta;
    
    // initialise grouping array
    group = vector< vector< vector<int> > >(n);
    for (int ind=0; ind<n; ind++) {
        group[ind] = vector< vector<int> >(loci);
        for (int l=0; l<loci; l++) {
            group[ind][l] = vector<int>(ploidy_vec[ind]);
        }
    }
    
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
    
    // initialise admix counts and frequencies
    admixCounts = vector< vector<int> >(n,vector<int>(K));
    admixCountsTotals = vector<int>(n);
    admixFreqs = vector< vector<double> >(n,vector<double>(K));
    
    // initialise objects for calculating assignment probabilities
    logProbVec = vector<double>(K);
    logProbVecSum = 0;  // (used in Qmatrix calculation)
    logProbVecMax = 0;
    probVec = vector<double>(K);
    probVecSum = 0;
    
    // initialise Qmatrices
    logQmatrix_gene_old = vector< vector< vector< vector<double> > > >(n);
    logQmatrix_gene_new = vector< vector< vector< vector<double> > > >(n);
    Qmatrix_gene_new = vector< vector< vector< vector<double> > > >(n);
    logQmatrix_gene_running = vector< vector< vector< vector<double> > > >(n);
    logQmatrix_gene = vector< vector< vector< vector<double> > > >(n);
    Qmatrix_gene = vector< vector< vector< vector<double> > > >(n);
    for (int ind=0; ind<n; ind++) {
        logQmatrix_gene_old[ind] = vector< vector< vector<double> > >(loci);
        logQmatrix_gene_new[ind] = vector< vector< vector<double> > >(loci);
        Qmatrix_gene_new[ind] = vector< vector< vector<double> > >(loci);
        logQmatrix_gene_running[ind] = vector< vector< vector<double> > >(loci);
        logQmatrix_gene[ind] = vector< vector< vector<double> > >(loci);
        Qmatrix_gene[ind] = vector< vector< vector<double> > >(loci);
        for (int l=0; l<loci; l++) {
            logQmatrix_gene_old[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K));
            logQmatrix_gene_new[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K));
            Qmatrix_gene_new[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K));
            logQmatrix_gene_running[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K));
            logQmatrix_gene[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K, log(double(0))));
            Qmatrix_gene[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K));
        }
    }
    
    // initialise objects for Hungarian algorithm
    costMat = vector< vector<double> >(K, vector<double>(K));
    bestPermOrder = vector<int>(K);
    
    edgesLeft = vector<int>(K);
    edgesRight = vector<int>(K);
    blockedLeft = vector<int>(K);
    blockedRight = vector<int>(K);
    
}

//------------------------------------------------
// particle_admixture::
// reset objects used in MCMC
void particle_admixture::reset(bool reset_Qmatrix_running) {
    
    // initialise group with random allocation
    vector<double> equalK(K,1/double(K));
    for (int ind=0; ind<n; ind++) {
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[l]; p++) {
                group[ind][l][p] = sample1(equalK,1.0);
            }
        }
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
                    thisGroup = group[ind][l][p];
                    alleleCounts[thisGroup-1][l][data[ind][l][p]-1]++;
                    alleleCountsTotals[thisGroup-1][l]++;
                }
            }
        }
    }
    
    // reset likelihoods
    logLikeGroup = 0;
    logLikeJoint = 0;
    
    // reset Qmatrices
    logQmatrix_gene_old = vector< vector< vector< vector<double> > > >(n);
    logQmatrix_gene_new = vector< vector< vector< vector<double> > > >(n);
    Qmatrix_gene_new = vector< vector< vector< vector<double> > > >(n);
    if (reset_Qmatrix_running) {
        logQmatrix_gene_running = vector< vector< vector< vector<double> > > >(n);
    }
    logQmatrix_gene = vector< vector< vector< vector<double> > > >(n);
    Qmatrix_gene = vector< vector< vector< vector<double> > > >(n);
    for (int ind=0; ind<n; ind++) {
        logQmatrix_gene_old[ind] = vector< vector< vector<double> > >(loci);
        logQmatrix_gene_new[ind] = vector< vector< vector<double> > >(loci);
        Qmatrix_gene_new[ind] = vector< vector< vector<double> > >(loci);
        if (reset_Qmatrix_running) {
            logQmatrix_gene_running[ind] = vector< vector< vector<double> > >(loci);
        }
        logQmatrix_gene[ind] = vector< vector< vector<double> > >(loci);
        Qmatrix_gene[ind] = vector< vector< vector<double> > >(loci);
        for (int l=0; l<loci; l++) {
            logQmatrix_gene_old[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K));
            logQmatrix_gene_new[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K));
            Qmatrix_gene_new[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K));
            if (reset_Qmatrix_running) {
                logQmatrix_gene_running[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K));
            }
            logQmatrix_gene[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K, log(double(0))));
            Qmatrix_gene[ind][l] = vector< vector<double> >(ploidy_vec[ind], vector<double>(K));
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
                thisData = data[ind][l][p];
        
                // subtract this gene copy from allele counts and admix counts
                if (thisData!=0) {   // if not missing data
                    thisGroup = group[ind][l][p];
                    alleleCounts[thisGroup-1][l][thisData-1]--;
                    alleleCountsTotals[thisGroup-1][l]--;
                    
                    admixCounts[ind][thisGroup-1]--;
                    admixCountsTotals[ind]--;
                }
                
                // calculate probability of this gene copy from all demes
                probVecSum = 0;
                for (int k=0; k<K; k++) {
                    if (thisData==0) {
                        probVec[k] = 1.0;
                    } else {
                        probVec[k] = double(alleleCounts[k][l][thisData-1]+lambda)/double(alleleCountsTotals[k][l]+J[l]*lambda);
                        if (beta!=1.0) {
                            probVec[k] = pow(probVec[k],beta);
                        }
                    }
                    probVec[k] *= double(admixCounts[ind][k]+alpha);  // (denominator of this expression is the same for all k, so is omitted)
                    probVecSum += probVec[k];
                }
                
                // resample grouping
                group[ind][l][p] = sample1(probVec, probVecSum);
                thisGroup = group[ind][l][p];
                
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
// choose best permutation of labels using method of Stephens (2000)
void particle_admixture::chooseBestLabelPermutation(globals &globals) {
    
    // calculate cost matrix from old and new Qmatrices
    for (int k1=0; k1<K; k1++) {
        for (int k2=0; k2<K; k2++) {
            costMat[k1][k2] = 0;
            for (int ind=0; ind<n; ind++) {
                for (int l=0; l<loci; l++) {
                    for (int p=0; p<ploidy_vec[ind]; p++) {
                        costMat[k1][k2] += Qmatrix_gene_new[ind][l][p][k1]*(logQmatrix_gene_new[ind][l][p][k1]-logQmatrix_gene_running[ind][l][p][k2]);
                    }
                }
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
        for (int ind=0; ind<n; ind++) {
            for (int l=0; l<loci; l++) {
                for (int p=0; p<ploidy_vec[ind]; p++) {
                    groupIndex++;
                    linearGroup[groupIndex] = bestPerm[linearGroup[groupIndex]-1]+1;
                    group[ind][l][p] = linearGroup[groupIndex];
                }
            }
        }
        
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
// particle_admixture::
// add logQmatrix_ind_new to logQmatrix_ind_running
void particle_admixture::updateQmatrix() {
    
    for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
            logQmatrix_ind_running[i][k] = logSum(logQmatrix_ind_running[i][k], logQmatrix_ind_new[i][k]);
        }
    }
}

//------------------------------------------------
// particle_admixture::
// store Qmatrix values
void particle_admixture::storeQmatrix() {
    
    // store individual-level Qmatrix
    for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
            logQmatrix_ind[i][k] = logSum(logQmatrix_ind[i][k], logQmatrix_ind_new[i][k]);
        }
    }
    
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
// draw allele frequencies given allele counts and lambda prior
void particle_admixture::drawFreqs() {
    
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
// particle_admixture::
// probability of data given grouping and allele frequencies
void particle_admixture::d_logLikeJoint() {
    
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
// particle_admixture::
// conditional probability of ith individual from kth deme (output in log space)
void particle_admixture::d_logLikeConditional(int i, int k) {
    
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



