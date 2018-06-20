//
//  MavericK
//  particle_noAdmixture.cpp
//
//  Created: Bob on 22/06/2016
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "particle_noAdmixture.h"

using namespace std;

extern vector< vector<double> > log_lookup;
extern vector<double> log_lookup_0;

//------------------------------------------------
// particle_noAdmixture::
// default constructor for class
particle_noAdmixture::particle_noAdmixture(){};

//------------------------------------------------
// particle_noAdmixture::
// informed constructor for class
particle_noAdmixture::particle_noAdmixture(globals &globals, int _K, double _beta_raised) {
    
    // copy some values over from globals object and arguments
    data = globals.data;
    K = _K;
    n = globals.n;
    loci = globals.loci;
    J = globals.J;
    ploidy_vec = globals.ploidy_vec;
    lambda = globals.lambda;
    beta_raised = _beta_raised;
    
    // initialise grouping vector
    group = vector<int>(n,1);
    group_propose = vector<int>(n);
    group_order = vector<int>(n);
    for (int i=0; i<n; i++) {
        group_order[i] = i+1;
    }
    
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
    
    // initialise objects for calculating assignment probabilities
    logProbVec = vector<double>(K);
    logProbVecMax = 0;
    probVec = vector<double>(K);
    probVecSum = 0;
    
    // initialise likelihoods
    logLikeGroup = 0;
    logLikeJoint = 0;
    
    // initialise Qmatrices
    logQmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    
}

//------------------------------------------------
// particle_noAdmixture::
// reset objects used in MCMC
void particle_noAdmixture::reset() {
    
    // initialise group with random allocation
    for (int i=0; i<n; i++) {
        group[i] = sample2(1,K);
    }
    
    // zero allele counts
    for (int k=0; k<K; k++) {
        for (int l=0; l<loci; l++) {
            alleleCountsTotals[k][l] = 0;
            for (int j=0; j<J[l]; j++) {
                alleleCounts[k][l][j] = 0;
            }
        }
    }
    
    // populate allele counts based on grouping
    for (int ind=0; ind<n; ind++) {
        for (int l=0; l<loci; l++) {
            int this_group = group[ind];
            for (int p=0; p<ploidy_vec[ind]; p++) {
                int this_data1 = data[ind][l][p];
                if (this_data1!=0) {
                    alleleCounts[this_group-1][l][this_data1-1]++;
                    alleleCountsTotals[this_group-1][l]++;
                }
            }
        }
    }
    
    // reset likelihoods
    logLikeGroup = 0;
    logLikeJoint = 0;
    
    // reset Qmatrices
    for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
            Qmatrix_ind[i][k] = 0;
            logQmatrix_ind[i][k] = 0;
        }
    }
    
}

//------------------------------------------------
// particle_noAdmixture::
// resample group allocation of all individuals by drawing from conditional posterior
void particle_noAdmixture::group_update() {
    
    // update group allocation for all individuals
    for (int ind=0; ind<n; ind++) {
        int this_group = group[ind];
        
        // subtract individual ind from allele counts
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                int this_data1 = data[ind][l][p];
                if (this_data1!=0) {   // if not missing data
                    alleleCounts[this_group-1][l][this_data1-1]--;
                    alleleCountsTotals[this_group-1][l]--;
                }
            }
        }
        
        // resample group allocation of individual ind
        if (beta_raised==0) {    // special case if beta==0 (draw from prior)
            this_group = sample2(1,K);
        } else {
            group_probs(ind);
            this_group = sample1(probVec, probVecSum);
        }
        group[ind] = this_group;
        
        // add individual ind to allele counts
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                int this_data1 = data[ind][l][p];
                if (this_data1!=0) {   // if not missing data
                    alleleCounts[this_group-1][l][this_data1-1]++;
                    alleleCountsTotals[this_group-1][l]++;
                }
            }
        }
        
    } // end loop over ind
    
}

//------------------------------------------------
// particle_noAdmixture::
// calculate conditional probability of individual ind coming from each deme
void particle_noAdmixture::group_probs(int ind) {
    
    for (int k=0; k<K; k++) {
        d_logLikeConditional(ind, k);   // recalculate logProbVec[k]
    }
    logProbVecMax = *max_element(begin(logProbVec),end(logProbVec));
    probVecSum = 0;
    for (int k=0; k<K; k++) {
        probVec[k] = exp(beta_raised*logProbVec[k] - beta_raised*logProbVecMax);
        probVecSum += probVec[k];
    }
    
}

//------------------------------------------------
// particle_noAdmixture::
// resample group allocation of all individuals by drawing from conditional posterior
void particle_noAdmixture::group_update_Klevel() {
    
    // if beta==0 then all moves accepted, so skip
    if (beta_raised==0) {
        return;
    }
    
    // choose two distinct random groups
    int K1 = sample2(1,K);
    int K2 = sample2(1,K);
    if (K2==K1) {
        K2++;
        if (K2>K) {
            K2 = 1;
        }
    }
    
    // search through observations in random order
    shuffle1(group_order);
    
    // subtract all individuals in these two groups
    double logLike_old = 0;
    double propose_logProb_old = 0;
    for (int i=(n-1); i>=0; i--) {  // loop backwards through group_order
        int ind = group_order[i]-1;
        int this_group = group[ind];
        if (this_group==K1 || this_group==K2) {
            
            // subtract individual ind from allele counts
            for (int l=0; l<loci; l++) {
                for (int p=0; p<ploidy_vec[ind]; p++) {
                    int this_data1 = data[ind][l][p];
                    if (this_data1!=0) {   // if not missing data
                        alleleCounts[this_group-1][l][this_data1-1]--;
                        alleleCountsTotals[this_group-1][l]--;
                    }
                }
            }
            
            // recalculate probability of allocation to K1 or K2
            d_logLikeConditional(ind, K1-1);   // recalculate logProbVec[K1]
            d_logLikeConditional(ind, K2-1);   // recalculate logProbVec[K2]
            logProbVecMax = (K1>K2) ? K1 : K2;
            probVec[K1-1] = exp(beta_raised*logProbVec[K1-1] - beta_raised*logProbVecMax);
            probVec[K2-1] = exp(beta_raised*logProbVec[K2-1] - beta_raised*logProbVecMax);
            probVecSum = probVec[K1-1] + probVec[K2-1];
            
            // calculate probability of chosen grouping
            propose_logProb_old += log(probVec[this_group-1]/probVecSum);
            logLike_old += beta_raised*logProbVec[this_group-1];
            
        }
    }   // end loop over individuals
    
    // propose a new grouping and calculate likelihood at same time
    double logLike_new = 0;
    double propose_logProb_new = 0;
    for (int i=0; i<n; i++) {  // loop forwards through group_order
        int ind = group_order[i]-1;
        int this_group = group[ind];
        if (this_group==K1 || this_group==K2) {
            
            // recalculate probability of allocation to K1 or K2
            d_logLikeConditional(ind, K1-1);   // recalculate logProbVec[K1]
            d_logLikeConditional(ind, K2-1);   // recalculate logProbVec[K2]
            logProbVecMax = (K1>K2) ? K1 : K2;
            probVec[K1-1] = exp(beta_raised*logProbVec[K1-1] - beta_raised*logProbVecMax);
            probVec[K2-1] = exp(beta_raised*logProbVec[K2-1] - beta_raised*logProbVecMax);
            probVecSum = probVec[K1-1] + probVec[K2-1];
            
            // resample group allocation of individual ind
            if (rbernoulli1(probVec[K1-1]/probVecSum)) {
                group_propose[ind] = K1;
            } else {
                group_propose[ind] = K2;
            }
            
            // calculate probability of new grouping
            propose_logProb_new += log(probVec[group_propose[ind]-1]/probVecSum);
            logLike_new += beta_raised*logProbVec[group_propose[ind]-1];
            
            // add individual ind to allele counts
            for (int l=0; l<loci; l++) {
                for (int p=0; p<ploidy_vec[ind]; p++) {
                    int this_data1 = data[ind][l][p];
                    if (this_data1!=0) {   // if not missing data
                        alleleCounts[group_propose[ind]-1][l][this_data1-1]++;
                        alleleCountsTotals[group_propose[ind]-1][l]++;
                    }
                }
            }
            
        }
    }   // end loop over individuals
    
    // Metropolis-Hastings step. If accept then stick with new grouping, otherwise revert back
    double MH_diff = (logLike_new - propose_logProb_new) - (logLike_old - propose_logProb_old);
    double rand1 = runif1(0,1);
    
    if (log(rand1)<MH_diff) {
        
        // accept move
        for (int ind=0; ind<n; ind++) {
            int this_group = group[ind];
            if (this_group==K1 || this_group==K2) {
                group[ind] = group_propose[ind];
            }
        }
        
    } else {
        
        // reject move
        for (int ind=0; ind<n; ind++) {
            int this_group = group[ind];
            if (this_group==K1 || this_group==K2) {
                
                for (int l=0; l<loci; l++) {
                    for (int p=0; p<ploidy_vec[ind]; p++) {
                        int this_data1 = data[ind][l][p];
                        if (this_data1!=0) {   // if not missing data
                            
                            // subtract new group
                            alleleCounts[group_propose[ind]-1][l][this_data1-1]--;
                            alleleCountsTotals[group_propose[ind]-1][l]--;
                            
                            // reinstate old group
                            alleleCounts[this_group-1][l][this_data1-1]++;
                            alleleCountsTotals[this_group-1][l]++;
                        }
                    }
                }
                
            }
        }   // loop over individuals
        
    }
    
}

//------------------------------------------------
// particle_noAdmixture::
// probability of data given grouping only, integrated over unknown allele frequencies
void particle_noAdmixture::d_logLikeGroup() {
    
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
// particle_noAdmixture::
// draw allele frequencies given allele counts and lambda prior
void particle_noAdmixture::drawFreqs() {
    
    // dirichlet distributed allele frequencies
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
// particle_noAdmixture::
// probability of data given grouping and allele frequencies
void particle_noAdmixture::d_logLikeJoint() {
    
    // calculate likelihood given known allele frequencies
    logLikeJoint = 0;
    double running = 1.0;
    for (int i=0; i<n; i++) {
        int this_group = group[i];
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[i]; p++) {
                int this_data1 = data[i][l][p];
                if (this_data1!=0) {
                    running *= alleleFreqs[this_group-1][l][this_data1-1];
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
// particle_noAdmixture::
// conditional probability of ith individual from kth deme (output in log space)
void particle_noAdmixture::d_logLikeConditional(int i, int k) {
    
    // reset logprob to zero
    logProbVec[k] = 0;
    
    // Likelihood calculation is highly tailored to different ploidy levels to increase speed and efficiency
    
    // IF DIPLOID
    if (ploidy_vec[i]==2) {
        
        // loop over loci
        for (unsigned int l=0; l<loci; l++) {
            
            // get both gene copies
            int this_data1 = data[i][l][0];
            int this_data2 = data[i][l][1];
            
            // if nethier gene copy missing then calculation is straightforward
            if (this_data1!=0 && this_data2!=0) {
                
                // get current alleleCounts and alleleCountsTotals
                int this_allele_counts = alleleCounts[k][l][this_data1-1];
                int this_allele_counts_totals = alleleCountsTotals[k][l];
                
                // if homozygote at this locus
                if (this_data1==this_data2) {
                    
                    // fast likelihood if within lookup table range
                    if ((this_allele_counts+1)<int(1e4) && (this_allele_counts_totals+1)<int(1e4)) {
                        logProbVec[k] += log_lookup_0[this_allele_counts] - log_lookup[this_allele_counts_totals][J[l]-1] + log_lookup_0[this_allele_counts+1] - log_lookup[this_allele_counts_totals+1][J[l]-1];
                    }
                    // otherwise slow likelihood
                    else {
                        logProbVec[k] += log((this_allele_counts + lambda)/double(this_allele_counts_totals + J[l]*lambda) * (this_allele_counts+1 + lambda)/double(this_allele_counts_totals+1 + J[l]*lambda));
                    }
                }
                // if heterozygote at this locus
                else {
                    
                    // fast likelihood if within lookup table range
                    if (this_allele_counts<int(1e4) && this_allele_counts_totals<int(1e4)) {
                        logProbVec[k] += log_lookup_0[this_allele_counts] - log_lookup[this_allele_counts_totals][J[l]-1];
                    }
                    // otherwise slow likelihood
                    else {
                        logProbVec[k] += log((this_allele_counts + lambda)/double(this_allele_counts_totals + J[l]*lambda));
                    }
                    
                    // get allele counts at second locus
                    this_allele_counts = alleleCounts[k][l][this_data2-1];
                    
                    // continue fast likelihood if within lookup table range
                    if (this_allele_counts<int(1e4) && (this_allele_counts_totals+1)<int(1e4)) {
                        logProbVec[k] += log_lookup_0[this_allele_counts] - log_lookup[this_allele_counts_totals+1][J[l]-1];
                    }
                    // otherwise slow likelihood
                    else {
                        logProbVec[k] += log((this_allele_counts + lambda)/double(this_allele_counts_totals+1 + J[l]*lambda));
                    }
                }
            }   // end if neither gene copy missing
            
            // if one or other gene copy missing
            else {
                
                // get current alleleCountsTotals
                int this_allele_counts_totals = alleleCountsTotals[k][l];
                
                // if first gene copy is present
                if (this_data1!=0) {
                    
                    // get current alleleCounts
                    int this_allele_counts = alleleCounts[k][l][this_data1-1];
                    
                    // fast likelihood if within lookup table range
                    if (this_allele_counts<int(1e4) && this_allele_counts_totals<int(1e4)) {
                        logProbVec[k] += log_lookup_0[this_allele_counts] - log_lookup[this_allele_counts_totals][J[l]-1];
                    }
                    // otherwise slow likelihood
                    else {
                        logProbVec[k] += log((this_allele_counts + lambda)/double(this_allele_counts_totals + J[l]*lambda));
                    }
                }
                // if second gene copy is present
                if (this_data2!=0) {
                    
                    // get current alleleCounts
                    int this_allele_counts = alleleCounts[k][l][this_data2-1];
                    
                    // fast likelihood if within lookup table range
                    if (this_allele_counts<int(1e4) && this_allele_counts_totals<int(1e4)) {
                        logProbVec[k] += log_lookup_0[this_allele_counts] - log_lookup[this_allele_counts_totals][J[l]-1];
                    }
                    // otherwise slow likelihood
                    else {
                        logProbVec[k] += log((this_allele_counts + lambda)/double(this_allele_counts_totals + J[l]*lambda));
                    }
                }
                
            }   // end if at least one gene copy missing
            
        } // end loop over loci
        
    } // end if diploid
    
    // IF HAPLOID
    else if (ploidy_vec[i]==1) {
        
        // loop over loci
        for (unsigned int l=0; l<loci; l++) {
            
            // get single gene copy
            int this_data1 = data[i][l][0];
            
            // if gene copy not missing
            if (this_data1!=0) {
                
                // get current alleleCounts and alleleCountsTotals
                int this_allele_counts = alleleCounts[k][l][this_data1-1];
                int this_allele_counts_totals = alleleCountsTotals[k][l];
                
                // fast likelihood if within lookup table range
                if (this_allele_counts<int(1e4) && this_allele_counts_totals<int(1e4)) {
                    logProbVec[k] += log_lookup_0[this_allele_counts] - log_lookup[this_allele_counts_totals][J[l]-1];
                }
                // otherwise slow likelihood
                else {
                    logProbVec[k] += log((this_allele_counts + lambda)/double(this_allele_counts_totals + J[l]*lambda));
                }
                
            }
            
        } // end loop over loci
        
    } // end if haploid
    
    // IF PLOIDY>2
    else {
        
        // loop over loci
        for (unsigned int l=0; l<loci; l++) {
            
            // get current alleleCountsTotals
            int this_allele_counts_totals = alleleCountsTotals[k][l];
            
            // loop through all gene copies
            for (unsigned int p=0; p<ploidy_vec[i]; p++) {
                
                // get this gene copy and check not missing
                int this_data1 = data[i][l][p];
                if (this_data1!=0) {
                    
                    // get current allele counts
                    int this_allele_counts = alleleCounts[k][l][this_data1-1];
                    
                    // fast likelihood if within lookup table range
                    if (this_allele_counts<int(1e4) && this_allele_counts_totals<int(1e4)) {
                        logProbVec[k] += log_lookup_0[this_allele_counts] - log_lookup[this_allele_counts_totals][J[l]-1];
                    }
                    // otherwise slow likelihood
                    else {
                        logProbVec[k] += log((this_allele_counts + lambda)/double(this_allele_counts_totals + J[l]*lambda));
                    }
                    
                    // increment allele counts
                    alleleCounts[k][l][this_data1-1] ++;
                    this_allele_counts_totals ++;
                }
            }
            // undo the temporary increment to allele counts used in likelihood calculation
            for (unsigned int p=0; p<ploidy_vec[i]; p++) {
                int this_data1 = data[i][l][p];
                if (this_data1!=0) {
                    alleleCounts[k][l][this_data1-1] --;
                }
            }
            
        } // end loop over loci
        
    } // end if ploidy>2
    
}



