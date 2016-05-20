//
//  MavericK
//  MCMCobject_noAdmixture.cpp
//
//  Created: Bob on 23/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "MCMCobject_noAdmixture.h"

using namespace std;

//------------------------------------------------
// MCMCobject_noAdmixture::
// constructor for class containing all elements required for MCMC under no-admixture model
MCMCobject_noAdmixture::MCMCobject_noAdmixture(globals &globals, int _Kindex, int _burnin, int _samples, int _thinning, double _beta) {
    
    // copy some values over from globals object
    outputQmatrix_pop_on = globals.outputQmatrix_pop_on;
    
    Kindex = _Kindex;
    K = globals.Kmin+Kindex;
    n = globals.n;
    loci = globals.loci;
    J = globals.J;
    ploidy_vec = globals.ploidy_vec;
    data = globals.data;
    lambda = globals.lambda;
    beta = _beta;
    uniquePops = globals.uniquePops;
    
    // create lookup table for log function
    int Jmax = *max_element(begin(J),end(J));
    log_lookup = vector< vector<double> >(int(1e4),vector<double>(Jmax+1));
    for (int i=0; i<int(1e4); i++) {
        for (int j=0; j<(Jmax+1); j++) {
            log_lookup[i][j] = log(double(i+j*lambda));
        }
    }
    
    burnin = _burnin;
    samples = _samples;
    thinning = _thinning;
    
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
    logQmatrix_ind_new = vector< vector<double> >(n, vector<double>(K));
    logQmatrix_ind_running = vector< vector<double> >(n, vector<double>(K));
    
    logQmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_pop = vector< vector<double> >(uniquePops.size(), vector<double>(K));
    
    // initialise objects for Hungarian algorithm
    costMat = vector< vector<double> >(K, vector<double>(K));
    bestPermOrder = vector<int>(K);
    
    edgesLeft = vector<int>(K);
    edgesRight = vector<int>(K);
    blockedLeft = vector<int>(K);
    blockedRight = vector<int>(K);
    
}

//------------------------------------------------
// MCMCobject_noAdmixture::
// reset objects used in MCMC
void MCMCobject_noAdmixture::reset(bool reset_Qmatrix_running) {
    
    // reset likelihoods
    logLikeGroup = 0;
    logLikeGroup_sum = 0;
    logLikeGroup_store = vector<double>(samples);
    logLikeGroup_sumSquared = 0;
    logLikeJoint = 0;
    logLikeJoint_sum = 0;
    logLikeJoint_sumSquared = 0;
    harmonic = log(double(0));
    
    // reset Qmatrices
    logQmatrix_ind_old = vector< vector<double> >(n, vector<double>(K));
    logQmatrix_ind_new = vector< vector<double> >(n, vector<double>(K));
    if (reset_Qmatrix_running) {
        logQmatrix_ind_running = vector< vector<double> >(n, vector<double>(K,-log(double(K))));
    }
    
    logQmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_pop = vector< vector<double> >(uniquePops.size(), vector<double>(K));
    
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
    
    // populate allele counts
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
    
}

//------------------------------------------------
// MCMCobject_noAdmixture::
// perform complete MCMC under no-admixture model
void MCMCobject_noAdmixture::perform_MCMC(globals &globals, bool drawAlleleFreqs, bool storeLoglike, bool fixLabels, bool outputLikelihood, bool outputPosteriorGrouping, int mainRep) {
    
    // perform MCMC
    int thinSwitch = 1;
    for (int rep=0; rep<(burnin+samples); rep++) {
        
        // thinning loop (becomes active after burn-in)
        for (int thin=0; thin<thinSwitch; thin++) {
            
            // update group allocation of all individuals
            group_update();
            
        }
        if (rep==burnin)
            thinSwitch = thinning;
        
        // if fix label-switching problem
        if (fixLabels) {
            // calculate logQmatrix_ind_new for this iteration
            produceQmatrix();
        
            // fix label-switching problem
            chooseBestLabelPermutation(globals, rep);
        
            // add logQmatrix_ind_new to logQmatrix_ind_running
            updateQmatrix(rep);
            
            // store Qmatrix values if no longer in burn-in
            if (rep>=burnin)
                storeQmatrix();
        }
        
        // calculate marginal likelihood
        d_logLikeGroup();
        
        // optionally draw allele frequencies and calculate joint likelihood
        if (drawAlleleFreqs) {
            drawFreqs();
            d_logLikeJoint();
        }
        
        // add likelihoods to running sums
        if (rep>=burnin) {
            logLikeGroup_sum += logLikeGroup;
            logLikeGroup_sumSquared += logLikeGroup*logLikeGroup;
            
            if (storeLoglike) {
                logLikeGroup_store[rep-burnin] = logLikeGroup;
            }

            harmonic = logSum(harmonic, -logLikeGroup);
            if (drawAlleleFreqs) {
                logLikeJoint_sum += logLikeJoint;
                logLikeJoint_sumSquared += logLikeJoint*logLikeJoint;
            }
        }
        
        // write to outputLikelihoods file
        if (outputLikelihood) {
            globals.outputLikelihood_fileStream << K << "," << mainRep+1 << "," << rep-burnin+1 << "," << logLikeGroup << "," << logLikeJoint << "\n";
            globals.outputLikelihood_fileStream.flush();
        }
        
        // write to outputPosteriorGrouping file
        if (outputPosteriorGrouping) {
            globals.outputPosteriorGrouping_fileStream << K << "," << mainRep+1 << "," << rep-burnin+1;
            for (int i=0; i<n; i++) {
                globals.outputPosteriorGrouping_fileStream << "," << group[i];
            }
            globals.outputPosteriorGrouping_fileStream << "\n";
            globals.outputPosteriorGrouping_fileStream.flush();
        }
            
    } // end of MCMC

    
    // finish off Qmatrices
    if (fixLabels) {
        
        // finish off individual level Qmatrix
        for (int i=0; i<n; i++) {
            for (int k=0; k<K; k++) {
                Qmatrix_ind[i][k] = exp(logQmatrix_ind[i][k] - log(double(samples)));
            }
        }
        
        // calculate population level Qmatrices
        if (outputQmatrix_pop_on) {
            for (int i=0; i<n; i++) {
                for (int k=0; k<K; k++) {
                    Qmatrix_pop[globals.pop_index[i]][k] += Qmatrix_ind[i][k];
                }
            }
            for (int i=0; i<int(uniquePops.size()); i++) {
                for (int k=0; k<K; k++) {
                    Qmatrix_pop[i][k] /= double(globals.uniquePop_counts[i]);
                }
            }
        }
    } // end if fixLabels
    
    // finish off harmonic mean
    harmonic = log(double(samples))-harmonic;
}

//------------------------------------------------
// MCMCobject_noAdmixture::
// resample group allocation of all individuals by drawing from conditional posterior
void MCMCobject_noAdmixture::group_update() {
    
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
// MCMCobject_noAdmixture::
// draw allele frequencies given allele counts and lambda prior
void MCMCobject_noAdmixture::drawFreqs() {

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
// MCMCobject_noAdmixture::
// choose best permutation of labels using method of Stephens (2000)
void MCMCobject_noAdmixture::chooseBestLabelPermutation(globals &globals, int rep) {
    
    // calculate cost matrix from old and new Qmatrices
    for (int k1=0; k1<K; k1++) {
        for (int k2=0; k2<K; k2++) {
            costMat[k1][k2] = 0;
            for (int i=0; i<n; i++) {
                costMat[k1][k2] += exp(logQmatrix_ind_new[i][k1])*(logQmatrix_ind_new[i][k1]-logQmatrix_ind_running[i][k2]);
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
// MCMCobject_noAdmixture::
// calculate logQmatrix_ind_new for this iteration
void MCMCobject_noAdmixture::produceQmatrix() {
    
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
    
}

//------------------------------------------------
// MCMCobject_noAdmixture::
// add logQmatrix_ind_new to logQmatrix_ind_running
void MCMCobject_noAdmixture::updateQmatrix(int &rep) {
    
    for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
            logQmatrix_ind_running[i][k] = logSum(logQmatrix_ind_running[i][k], logQmatrix_ind_new[i][k]);
        }
    }
}

//------------------------------------------------
// MCMCobject_noAdmixture::
// store Qmatrix values
void MCMCobject_noAdmixture::storeQmatrix() {
    
    // store individual-level Qmatrix
    for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
            logQmatrix_ind[i][k] = logSum(logQmatrix_ind[i][k], logQmatrix_ind_new[i][k]);
        }
    }
    
}

//------------------------------------------------
// MCMCobject_noAdmixture::
// conditional probability of ith individual from kth deme (output in log space)
void MCMCobject_noAdmixture::d_logLikeConditional(int i, int k) {
    
    // calculate conditional probability of data
    logProbVec[k] = 0;
    int d, a, a_t;  // for making temporary copies of data, alleleCounts, and alleleCountsTotals respectively
    for (unsigned int l=0; l<loci; l++) {
        a_t = alleleCountsTotals[k][l];
        for (unsigned int p=0; p<ploidy_vec[i]; p++) {
            d = data[i][l][p];
            a = alleleCounts[k][l][d-1];
            if (d!=0) {
                if ((a<int(1e4)) && (a_t<int(1e4))) {
                    logProbVec[k] += log_lookup[a][1]-log_lookup[a_t][J[l]];
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

//------------------------------------------------
// MCMCobject_noAdmixture::
// probability of data given grouping only, integrated over unknown allele frequencies
void MCMCobject_noAdmixture::d_logLikeGroup() {
    
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
// MCMCobject_noAdmixture::
// probability of data given grouping and allele frequencies
void MCMCobject_noAdmixture::d_logLikeJoint() {
    
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
