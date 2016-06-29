//
//  MavericK
//  EM_algorithm.cpp
//
//  Created: Bob on 25/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "EM_algorithm.h"

using namespace std;

//------------------------------------------------
// EM algorithm under no-admixture model
void EM_noAdmix(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // create empty objects
    vector< vector<double> > logProbMat(globals.n,vector<double>(K));
    double logProbSum;
    vector< vector< vector<double> > > logAlleleFreqs(K);
    double logAlleleFreqsSum;
    double logLike=0;
    double logLike_this_i=0;
    double logLike_this_ik=0;
    
    // initialise logAlleleFreqs and zero global maximum likelihood allele frequencies
    globals.max_alleleFreqs = vector< vector< vector<double> > >(K);
    for (int k=0; k<K; k++) {
        logAlleleFreqs[k] = vector< vector<double> >(globals.loci);
        globals.max_alleleFreqs[k] = vector< vector<double> >(globals.loci);
        for (int l=0; l<globals.loci; l++) {
            logAlleleFreqs[k][l] = vector<double>(globals.J[l]);
            globals.max_alleleFreqs[k][l] = vector<double>(globals.J[l]);
        }
    }
    
    // repeat EM multiple times
    globals.maxLike = log(double(0));
    for (int EMrep=0; EMrep<globals.EMrepeats; EMrep++) {
        
        // initialise random allele frequencies
        for (int k=0; k<K; k++) {
            for (int l=0; l<globals.loci; l++) {
                logAlleleFreqsSum = log(double(0));
                for (int j=0; j<globals.J[l]; j++) {
                    logAlleleFreqs[k][l][j] = log(runif1(0.1,0.9));
                    logAlleleFreqsSum = logSum(logAlleleFreqsSum, logAlleleFreqs[k][l][j]);
                }
                for (int j=0; j<globals.J[l]; j++) {
                    logAlleleFreqs[k][l][j] -= logAlleleFreqsSum;
                }
            }
        }
        
        // EM iterations
        for (int EMiter=0; EMiter<globals.EMiterations; EMiter++) {
            
            // calculate assignment probability matrix
            for (int i=0; i<globals.n; i++) {
                logProbSum = log(double(0));
                for (int k=0; k<K; k++) {
                    logProbMat[i][k] = 0;
                    for (int l=0; l<globals.loci; l++) {
                        for (int p=0; p<globals.ploidy_vec[i]; p++) {
                            if (globals.data[i][l][p]!=0) {
                                logProbMat[i][k] += logAlleleFreqs[k][l][globals.data[i][l][p]-1];
                            }
                        }
                    }
                    logProbSum = logSum(logProbSum, logProbMat[i][k]);
                }
                for (int k=0; k<K; k++) {
                    logProbMat[i][k] -= logProbSum;
                }
            }
            
            // zero allele frequencies
            for (int k=0; k<K; k++) {
                for (int l=0; l<globals.loci; l++) {
                    for (int j=0; j<globals.J[l]; j++) {
                        logAlleleFreqs[k][l][j] = log(double(0));
                    }
                }
            }
            
            // add probability-weighted contribution to allele frequencies
            for (int i=0; i<globals.n; i++) {
                for (int l=0; l<globals.loci; l++) {
                    for (int p=0; p<globals.ploidy_vec[i]; p++) {
                        if (globals.data[i][l][p]!=0) {
                            for (int k=0; k<K; k++) {
                                logAlleleFreqs[k][l][globals.data[i][l][p]-1] = logSum(logAlleleFreqs[k][l][globals.data[i][l][p]-1], logProbMat[i][k]);
                            }
                        }
                    }
                }
            }
            
            // normalise allele frequencies
            for (int k=0; k<K; k++) {
                for (int l=0; l<globals.loci; l++) {
                    logAlleleFreqsSum = log(double(0));
                    for (int j=0; j<globals.J[l]; j++) {
                        // control for underflow
                        if (logAlleleFreqs[k][l][j] < -1e300)
                            logAlleleFreqs[k][l][j] = -1e300;
                        logAlleleFreqsSum = logSum(logAlleleFreqsSum, logAlleleFreqs[k][l][j]);
                    }
                    for (int j=0; j<globals.J[l]; j++) {
                        logAlleleFreqs[k][l][j] -= logAlleleFreqsSum;
                    }
                }
            }
            
        } // end of EM iterations loop
        
        // calculate log-likelihood
        logLike = 0;
        for (int i=0; i<globals.n; i++) {
            logLike_this_i = log(double(0));
            for (int k=0; k<K; k++) {
                logLike_this_ik = 0;
                for (int l=0; l<globals.loci; l++) {
                    for (int p=0; p<globals.ploidy_vec[i]; p++) {
                        if (globals.data[i][l][p]!=0) {
                            logLike_this_ik += logAlleleFreqs[k][l][globals.data[i][l][p]-1];
                        }
                    }
                }
                logLike_this_i = logSum(logLike_this_i, logLike_this_ik - log(double(K)));
            }
            logLike += logLike_this_i;
        }
        
        // replace current values if more likely
        if (logLike>globals.maxLike) {
            globals.maxLike = logLike;
            for (int k=0; k<K; k++) {
                for (int l=0; l<globals.loci; l++) {
                    for (int j=0; j<globals.J[l]; j++) {
                        globals.max_alleleFreqs[k][l][j] = exp(logAlleleFreqs[k][l][j]);
                    }
                }
            }
        }
        
    } // end of repeat EM loop
    
    // calculate model comparison statistics
    if (globals.outputComparisonStatistics_on) {
        
        int freeParameters = K*(sum(globals.J)-globals.loci);
        globals.AIC[Kindex] = 2*freeParameters - 2*globals.maxLike;
        globals.BIC[Kindex] = freeParameters*log(double(globals.n)) - 2*globals.maxLike;
        globals.DIC_Spiegelhalter[Kindex] = -4*globals.structure_loglike_mean[Kindex] + 2*globals.maxLike;
        globals.DIC_Gelman[Kindex] = -2*globals.structure_loglike_mean[Kindex] + 2*globals.structure_loglike_var[Kindex];
        
    }
    
}

//------------------------------------------------
// EM algorithm under admixture model
void EM_admix(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // create empty objects
    int geneCopies = sum(globals.ploidy_vec)*globals.loci;
    vector< vector<double> > logProbMat(geneCopies,vector<double>(K,-log(double(K))));
    double logProbSum;
    vector< vector< vector<double> > > logAlleleFreqs(K);
    double logAlleleFreqsSum;
    vector< vector<double> > logAdmixFreqs(globals.n,vector<double>(K));
    double logLike=0;
    double logLike_this_i=0;
    double logLike_this_ik=0;
    
    // initialise logAlleleFreqs and logAdmixFreqs and zero global maximum likelihood allele frequencies
    globals.max_alleleFreqs = vector< vector< vector<double> > >(K);
    globals.max_admixFreqs = vector< vector<double> >(globals.n,vector<double>(K));
    for (int k=0; k<K; k++) {
        logAlleleFreqs[k] = vector< vector<double> >(globals.loci);
        globals.max_alleleFreqs[k] = vector< vector<double> >(globals.loci);
        for (int l=0; l<globals.loci; l++) {
            logAlleleFreqs[k][l] = vector<double>(globals.J[l]);
            globals.max_alleleFreqs[k][l] = vector<double>(globals.J[l]);
        }
    }
    
    // repeat EM multiple times
    globals.maxLike = log(double(0));
    for (int EMrep=0; EMrep<globals.EMrepeats; EMrep++) {
        
        // initialise random allele frequencies
        for (int k=0; k<K; k++) {
            for (int l=0; l<globals.loci; l++) {
                logAlleleFreqsSum = log(double(0));
                for (int j=0; j<globals.J[l]; j++) {
                    logAlleleFreqs[k][l][j] = log(runif1(0.1,0.9));
                    logAlleleFreqsSum = logSum(logAlleleFreqsSum, logAlleleFreqs[k][l][j]);
                }
                for (int j=0; j<globals.J[l]; j++) {
                    logAlleleFreqs[k][l][j] -= logAlleleFreqsSum;
                }
            }
        }
        
        // EM iterations
        for (int EMiter=0; EMiter<globals.EMiterations; EMiter++) {
            
            // loop through all individuals
            int groupIndex = -1;
            for (int i=0; i<globals.n; i++) {
                
                // zero admixture freqs for this individual
                for (int k=0; k<K; k++) {
                    logAdmixFreqs[i][k] = log(double(0));
                }
                
                // calculate probability matrix and contribute to ML admixture freqs
                for (int l=0; l<globals.loci; l++) {
                    for (int p=0; p<globals.ploidy_vec[i]; p++) {
                        groupIndex++;
                        if (globals.data[i][l][p]!=0) {
                            logProbSum = log(double(0));
                            for (int k=0; k<K; k++) {
                                logProbMat[groupIndex][k] = logAlleleFreqs[k][l][globals.data[i][l][p]-1];
                                logProbSum = logSum(logProbSum,logProbMat[groupIndex][k]);
                            }
                            for (int k=0; k<K; k++) {
                                logProbMat[groupIndex][k] -= logProbSum;
                            }
                        }
                        for (int k=0; k<K; k++) {
                            logAdmixFreqs[i][k] = logSum(logAdmixFreqs[i][k],logProbMat[groupIndex][k]-log(double(globals.loci))-log(double(globals.ploidy_vec[i])));
                        }
                    }
                }
            } // i loop
            
            // zero allele frequencies
            for (int k=0; k<K; k++) {
                for (int l=0; l<globals.loci; l++) {
                    for (int j=0; j<globals.J[l]; j++) {
                        logAlleleFreqs[k][l][j] = log(double(0));
                    }
                }
            }
            
            // add contribution to allele frequencies
            groupIndex = -1;
            for (int i=0; i<globals.n; i++) {
                for (int l=0; l<globals.loci; l++) {
                    for (int p=0; p<globals.ploidy_vec[i]; p++) {
                        groupIndex++;
                        if (globals.data[i][l][p]!=0) {
                            for (int k=0; k<K; k++) {
                                logAlleleFreqs[k][l][globals.data[i][l][p]-1] = logSum(logAlleleFreqs[k][l][globals.data[i][l][p]-1], logAdmixFreqs[i][k]+logProbMat[groupIndex][k]);
                            }
                        }
                    }
                }
            }
            
            // normalise allele frequencies
            for (int k=0; k<K; k++) {
                for (int l=0; l<globals.loci; l++) {
                    logProbSum = log(double(0));
                    for (int j=0; j<globals.J[l]; j++) {
                        logProbSum = logSum(logProbSum, logAlleleFreqs[k][l][j]);
                    }
                    for (int j=0; j<globals.J[l]; j++) {
                        logAlleleFreqs[k][l][j] -= logProbSum;
                    }
                }
            }
            
        } // end of EM iterations loop
        
        
        logLike = 0;
        for (int i=0; i<globals.n; i++) {
            logLike_this_i = log(double(0));
            for (int k=0; k<K; k++) {
                logLike_this_ik = 0;
                for (int l=0; l<globals.loci; l++) {
                    for (int p=0; p<globals.ploidy_vec[i]; p++) {
                        if (globals.data[i][l][p]!=0) {
                            logLike_this_ik += logAlleleFreqs[k][l][globals.data[i][l][p]-1];
                        }
                    }
                }
                logLike_this_i = logSum(logLike_this_i, logLike_this_ik - log(double(K)));
            }
            logLike += logLike_this_i;
        }
        
        // calculate log-likelihood
        double logLike = 0;
        double temp;
        for (int i=0; i<globals.n; i++) {
            for (int l=0; l<globals.loci; l++) {
                for (int p=0; p<globals.ploidy; p++) {
                    if (globals.data[i][l][p]!=0) {
                        temp = log(double(0));
                        for (int k=0; k<K; k++) {
                            temp = logSum(temp, logAdmixFreqs[i][k]+logAlleleFreqs[k][l][globals.data[i][l][p]-1]);
                        }
                        logLike += temp;
                    }
                }
            }
        }
        
        // replace current values if more likely
        if (logLike>globals.maxLike) {
            globals.maxLike = logLike;
            for (int k=0; k<K; k++) {
                for (int l=0; l<globals.loci; l++) {
                    for (int j=0; j<globals.J[l]; j++) {
                        globals.max_alleleFreqs[k][l][j] = exp(logAlleleFreqs[k][l][j]);
                    }
                }
            }
            for (int i=0; i<globals.n; i++) {
                for (int k=0; k<K; k++) {
                    globals.max_admixFreqs[i][k] = exp(logAdmixFreqs[i][k]);
                }
            }
        }
        
    } // end of repeat EM loop
    
    // calculate model comparison statistics
    if (globals.outputComparisonStatistics_on) {
        
        int freeParameters = K*(sum(globals.J)-globals.loci) + globals.n*(K-1);
        globals.AIC[Kindex] = 2*freeParameters - 2*globals.maxLike;
        globals.BIC[Kindex] = freeParameters*log(double(globals.geneCopies)) - 2*globals.maxLike;
        globals.DIC_Spiegelhalter[Kindex] = -4*globals.structure_loglike_mean[Kindex] + 2*globals.maxLike;
        globals.DIC_Gelman[Kindex] = -2*globals.structure_loglike_mean[Kindex] + 2*globals.structure_loglike_var[Kindex];
        
    }
    
}