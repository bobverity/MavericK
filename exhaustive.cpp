//
//  MavericK
//  exhaustive.cpp
//
//  Created: Bob on 23/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "exhaustive.h"

using namespace std;

//------------------------------------------------
// exhaustive analysis under no-admixture model
void exhaustive_noAdmix(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // initialise group
    vector<int> group(globals.n,1);
    vector<int> newGroup = group;
    
    // initialise alleleCounts array
    vector< vector< vector<int> > > alleleCounts(K);
    vector< vector<int> > alleleCountsTotals(K);
    for (int k=0; k<K; k++) {
        alleleCounts[k] = vector< vector<int> >(globals.loci);
        alleleCountsTotals[k] = vector<int>(globals.loci);
        for (int l=0; l<globals.loci; l++) {
            alleleCounts[k][l] = vector<int>(globals.J[l]);
        }
    }
    
    // populate alleleCounts array and calculate starting likelihood
    double logLike=0;
    for (int ind=0; ind<globals.n; ind++) {
        for (int l=0; l<globals.loci; l++) {
            for (int p=0; p<globals.ploidy_vec[ind]; p++) {
                if (globals.data[ind][l][p]!=0) {
                    logLike += log(alleleCounts[group[ind]-1][l][globals.data[ind][l][p]-1] + globals.lambda) - log(alleleCountsTotals[group[ind]-1][l] + globals.J[l]*globals.lambda);
                    alleleCounts[group[ind]-1][l][globals.data[ind][l][p]-1]++;
                    alleleCountsTotals[group[ind]-1][l]++;
                }
            }
        }
    }
    double combined_logLike = logLike -globals.n*log(double(K));
    combined_logLike += log(double(K));
    
    // loop through all possible groupings
    int uniques = 0;  // number of unique elements in group. Forms part of combinatorial constant needed to calculate correct likelihood when only searching over unique partitions (see Notes.c for details)
    double logLike_total = combined_logLike;
    while (group[0]!=0) {
        
        increment_restrictedGrowth(newGroup,K);
        if (newGroup[0]==0)
            break;
        
        // loop through all individuals
        for (int ind=0; ind<globals.n; ind++) {
            
            // if group allocation for this individual has changed
            if (newGroup[ind]!=group[ind]) {
                
                // subtract old group from allele counts and likelihood
                for (int l=0; l<globals.loci; l++) {
                    for (int p=0; p<globals.ploidy_vec[ind]; p++) {
                        if (globals.data[ind][l][p]!=0) {
                            alleleCountsTotals[group[ind]-1][l]--;
                            alleleCounts[group[ind]-1][l][globals.data[ind][l][p]-1]--;
                            logLike -= log(alleleCounts[group[ind]-1][l][globals.data[ind][l][p]-1] + globals.lambda) - log(alleleCountsTotals[group[ind]-1][l] + globals.J[l]*globals.lambda);
                        }
                    }
                }
                
                // add new group to allele counts and likelihood
                for (int l=0; l<globals.loci; l++) {
                    for (int p=0; p<globals.ploidy_vec[ind]; p++) {
                        if (globals.data[ind][l][p]!=0) {
                            logLike += log(alleleCounts[newGroup[ind]-1][l][globals.data[ind][l][p]-1] + globals.lambda) - log(alleleCountsTotals[newGroup[ind]-1][l] + globals.J[l]*globals.lambda);
                            alleleCounts[newGroup[ind]-1][l][globals.data[ind][l][p]-1]++;
                            alleleCountsTotals[newGroup[ind]-1][l]++;
                        }
                    }
                }
                
            }
        }
        
        // update group
        group = newGroup;
        
        // some combinatorics to account for searching unique partitions only (see Notes.c for details)
        combined_logLike = logLike -globals.n*log(double(K));
        uniques = *max_element(begin(group),end(group));
        combined_logLike += my_lgamma(double(K)+1)-my_lgamma(double(K)-double(uniques)+1);
        
        // add to running sum
        logLike_total = logSum(logLike_total, combined_logLike);
    
    }
    
    globals.logEvidence_exhaustive[Kindex] = logLike_total;
    
}

//------------------------------------------------
// exhaustive analysis under admixture model (alpha fixed or variable)
void exhaustive_admix(globals &globals, int Kindex) {
    
    // alpha fixed
    if (globals.fixAlpha_on) {
        
        globals.logEvidence_exhaustive[Kindex] = exhaustive_admix_fixedAlpha(globals, Kindex, globals.alpha[Kindex]);
        
    // alpha variable (integrate over alpha by brute force)
    } else {
        
        double logEvidence_total = log(double(0));
        for (int i=0; i<100; i++) {
            double alpha = (i+1)/10.0;
            logEvidence_total = logSum(logEvidence_total, exhaustive_admix_fixedAlpha(globals, Kindex, alpha));
        }
        logEvidence_total -= log(100.0);
        globals.logEvidence_exhaustive[Kindex] = logEvidence_total;
        
    }
    
}

//------------------------------------------------
// exhaustive analysis under admixture model for given alpha
double exhaustive_admix_fixedAlpha(globals &globals, int Kindex, double alpha) {
    int K = globals.Kmin+Kindex;
    
    // initialise group
    int geneCopies = sum(globals.ploidy_vec)*globals.loci;
    vector<int> group(geneCopies,1);
    vector<int> newGroup = group;
    
    // initialise alleleCounts array
    vector< vector< vector<int> > > alleleCounts(K);
    vector< vector<int> > alleleCountsTotals(K);
    for (int k=0; k<K; k++) {
        alleleCounts[k] = vector< vector<int> >(globals.loci);
        alleleCountsTotals[k] = vector<int>(globals.loci);
        for (int l=0; l<globals.loci; l++) {
            alleleCounts[k][l] = vector<int>(globals.J[l]);
        }
    }
    
    // initialise admixCounts array
    vector< vector<int> > admixCounts(globals.n, vector<int>(K));
    vector<int> admixCountsTotals(globals.n);
    
    // populate alleleCounts and admixCounts arrays and calculate starting likelihood
    double logLike=0;
    int linearIndex=-1;
    for (int ind=0; ind<globals.n; ind++) {
        for (int l=0; l<globals.loci; l++) {
            for (int p=0; p<globals.ploidy_vec[ind]; p++) {
                linearIndex++;
                if (globals.data[ind][l][p]!=0) {
                    logLike += log(alleleCounts[group[linearIndex]-1][l][globals.data[ind][l][p]-1] + globals.lambda) - log(alleleCountsTotals[group[linearIndex]-1][l] + globals.J[l]*globals.lambda);
                    alleleCounts[group[linearIndex]-1][l][globals.data[ind][l][p]-1]++;
                    alleleCountsTotals[group[linearIndex]-1][l]++;
                    
                    logLike += log(admixCounts[ind][group[linearIndex]-1] + alpha) - log(admixCountsTotals[ind] + K*alpha);
                    admixCounts[ind][group[linearIndex]-1]++;
                    admixCountsTotals[ind]++;
                }
            }
        }
    }
    double combined_logLike = logLike + log(double(K));
    
    // loop through all possible groupings
    int uniques = 0;  // number of unique elements in group. Forms part of combinatorial constant needed to calculate correct likelihood when only searching over unique partitions (see Notes.c for details)
    double logLike_total = combined_logLike;
    while (group[0]!=0) {
        
        increment_restrictedGrowth(newGroup,K);
        if (newGroup[0]==0)
            break;
        
        // loop through all gene copies
        linearIndex=-1;
        for (int ind=0; ind<globals.n; ind++) {
            for (int l=0; l<globals.loci; l++) {
                for (int p=0; p<globals.ploidy_vec[ind]; p++) {
                    linearIndex++;
                    
                    // if group allocation for this gene copy has changed
                    if (newGroup[linearIndex]!=group[linearIndex]) {
                        
                        // subtract old group from allele counts, admix counts and likelihood
                        if (globals.data[ind][l][p]!=0) {
                            alleleCountsTotals[group[linearIndex]-1][l]--;
                            alleleCounts[group[linearIndex]-1][l][globals.data[ind][l][p]-1]--;
                            logLike -= log(alleleCounts[group[linearIndex]-1][l][globals.data[ind][l][p]-1] + globals.lambda) - log(alleleCountsTotals[group[linearIndex]-1][l] + globals.J[l]*globals.lambda);
                            
                            admixCountsTotals[ind]--;
                            admixCounts[ind][group[linearIndex]-1]--;
                            logLike -= log(admixCounts[ind][group[linearIndex]-1] + alpha) - log(admixCountsTotals[ind] + K*alpha);
                        }
                        
                        // add new group to allele counts, admix counts and likelihood
                        if (globals.data[ind][l][p]!=0) {
                            logLike += log(alleleCounts[newGroup[linearIndex]-1][l][globals.data[ind][l][p]-1] + globals.lambda) - log(alleleCountsTotals[newGroup[linearIndex]-1][l] + globals.J[l]*globals.lambda);
                            alleleCounts[newGroup[linearIndex]-1][l][globals.data[ind][l][p]-1]++;
                            alleleCountsTotals[newGroup[linearIndex]-1][l]++;
                            
                            logLike += log(admixCounts[ind][newGroup[linearIndex]-1] + alpha) - log(admixCountsTotals[ind] + K*alpha);
                            admixCounts[ind][newGroup[linearIndex]-1]++;
                            admixCountsTotals[ind]++;
                        }
                        
                    }
                    
                } // p
            } // l
        } // ind
        
        // update group
        group = newGroup;
        
        // some combinatorics to account for searching unique partitions only (see Notes.c for details)
        uniques = *max_element(begin(group),end(group));
        combined_logLike = logLike + my_lgamma(double(K)+1)-my_lgamma(double(K)-double(uniques)+1);
        
        // add to running sum
        logLike_total = logSum(logLike_total, combined_logLike);
    }
    
    return(logLike_total);
}