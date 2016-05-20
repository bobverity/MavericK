//
//  MavericK
//  mainMCMC.cpp
//
//  Created: Bob on 23/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "mainMCMC.h"

using namespace std;

//------------------------------------------------
// main Structure MCMC under no-admixture model, repeated multiple times
void mainMCMC_noAdmixture(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // define temporary Qmatrix objects for calculating mean and sd over mainRepeats
    vector< vector< vector<double> > > Qmatrix_ind_allReps(globals.n, vector< vector<double> >(K, vector<double>(globals.mainRepeats)));
    vector< vector< vector<double> > > Qmatrix_pop_allReps(globals.uniquePops.size(), vector< vector<double> >(K, vector<double>(globals.mainRepeats)));
    
    // define MCMC object
    MCMCobject_noAdmixture mainMCMC(globals, Kindex, globals.mainBurnin, globals.mainSamples, globals.mainThinning, 1.0);
    
    // repeat analysis multiple times
    for (int mainRep=0; mainRep<globals.mainRepeats; mainRep++) {
        coutAndLog("  analysis "+to_string((long long)mainRep+1)+" of "+to_string((long long)globals.mainRepeats)+"\n", globals.outputLog_on, globals.outputLog_fileStream);
        
        // perform MCMC
        if (mainRep==0) {
            mainMCMC.reset(true);
        } else {
            mainMCMC.reset(false);
        }
        mainMCMC.perform_MCMC(globals, true, false, globals.fixLabels_on, globals.outputLikelihood_on, globals.outputPosteriorGrouping_on, mainRep);
        
        // save Qmatrix values
        if (globals.fixLabels_on) {
            for (int i=0; i<globals.n; i++) {
                for (int k=0; k<K; k++) {
                    Qmatrix_ind_allReps[i][k][mainRep] = mainMCMC.Qmatrix_ind[i][k];
                }
            }
            if (globals.outputQmatrix_pop_on) {
                for (int i=0; i<int(globals.uniquePops.size()); i++) {
                    for (int k=0; k<K; k++) {
                        Qmatrix_pop_allReps[i][k][mainRep] = mainMCMC.Qmatrix_pop[i][k];
                    }
                }
            }
        }
        
        // save harmonic mean of marginal likelihoods
        globals.logEvidence_harmonic[Kindex][mainRep] = mainMCMC.harmonic;
        
        // calculate Structure estimator from joint likelihoods
        globals.structure_loglike_mean[Kindex][mainRep] = mainMCMC.logLikeJoint_sum/double(globals.mainSamples);
        globals.structure_loglike_var[Kindex][mainRep] = mainMCMC.logLikeJoint_sumSquared/double(globals.mainSamples) - globals.structure_loglike_mean[Kindex][mainRep]*globals.structure_loglike_mean[Kindex][mainRep];
        globals.logEvidence_structure[Kindex][mainRep] = globals.structure_loglike_mean[Kindex][mainRep] - 0.5*globals.structure_loglike_var[Kindex][mainRep];
        
        // calculate elements of Evanno's delta K for the *previous* iteration
        if (Kindex>0)
            globals.L_1[Kindex][mainRep] = globals.logEvidence_structure[Kindex][mainRep]-globals.logEvidence_structure[Kindex-1][mainRep];
        if (Kindex>1) {
            globals.L_2[Kindex-1][mainRep] = abs(globals.L_1[Kindex][mainRep]-globals.L_1[Kindex-1][mainRep]);
        }
        
    }
    
    // save Qmatrices
    if (globals.fixLabels_on) {
        for (int i=0; i<globals.n; i++) {
            for (int k=0; k<K; k++) {
                globals.Qmatrix_ind[Kindex][i][k] = mean(Qmatrix_ind_allReps[i][k]);
                globals.QmatrixError_ind[Kindex][i][k] = sqrt(var(Qmatrix_ind_allReps[i][k])/double(globals.mainRepeats));
            }
        }
        if (globals.outputQmatrix_pop_on) {
            for (int i=0; i<int(globals.uniquePops.size()); i++) {
                for (int k=0; k<K; k++) {
                    globals.Qmatrix_pop[Kindex][i][k] = mean(Qmatrix_pop_allReps[i][k]);
                    globals.QmatrixError_pop[Kindex][i][k] = sqrt(var(Qmatrix_pop_allReps[i][k]));
                }
            }
        }
    }
    
    // calculate grand mean and standard error of harmonic mean estimator
    globals.logEvidence_harmonic_grandMean[Kindex] = mean(globals.logEvidence_harmonic[Kindex]);
    if (globals.mainRepeats>1) {
        globals.logEvidence_harmonic_grandSE[Kindex] = sqrt(var(globals.logEvidence_harmonic[Kindex], 0, true)/globals.mainRepeats);
    }
    
    // calculate grand mean and standard error of Structure estimator
    globals.logEvidence_structure_grandMean[Kindex] = mean(globals.logEvidence_structure[Kindex]);
    if (globals.mainRepeats>1) {
        globals.logEvidence_structure_grandSE[Kindex] = sqrt(var(globals.logEvidence_structure[Kindex], 0, true)/globals.mainRepeats);
    }
    
    // finalise calculation of Evanno's delta K for *previous* iteration
    if (Kindex>1) {
        globals.delta_K[Kindex-1] = mean(globals.L_2[Kindex-1])/sqrt(var(globals.logEvidence_structure[Kindex-1], 0, true));
    }
    
}

//------------------------------------------------
// main Structure MCMC under admixture model, repeated multiple times
void mainMCMC_admixture(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // define temporary Qmatrix objects for calculating mean and sd over mainRepeats
    vector< vector< vector<double> > > Qmatrix_gene_allReps(globals.geneCopies, vector< vector<double> >(K, vector<double>(globals.mainRepeats)));
    vector< vector< vector<double> > > Qmatrix_ind_allReps(globals.n, vector< vector<double> >(K, vector<double>(globals.mainRepeats)));
    vector< vector< vector<double> > > Qmatrix_pop_allReps(globals.uniquePops.size(), vector< vector<double> >(K, vector<double>(globals.mainRepeats)));
    
    // define MCMC object
    MCMCobject_admixture mainMCMC(globals, Kindex, globals.mainBurnin, globals.mainSamples, globals.mainThinning, 1.0);
    
    // repeat analysis multiple times
    for (int mainRep=0; mainRep<globals.mainRepeats; mainRep++) {
        coutAndLog("  analysis "+to_string((long long)mainRep+1)+" of "+to_string((long long)globals.mainRepeats)+"\n", globals.outputLog_on, globals.outputLog_fileStream);
        
        // perform MCMC
        if (mainRep==0) {
            mainMCMC.reset(true);
        } else {
            mainMCMC.reset(false);
        }
		mainMCMC.perform_MCMC(globals, true, false, globals.fixLabels_on, globals.outputLikelihood_on, globals.outputPosteriorGrouping_on, mainRep);
        
        // save Qmatrix values
        if (globals.fixLabels_on) {
            for (int i=0; i<globals.geneCopies; i++) {
                for (int k=0; k<K; k++) {
                    Qmatrix_gene_allReps[i][k][mainRep] = mainMCMC.Qmatrix_gene[i][k];
                    /*
                    if (k!=0)
                        globals.junk_fileStream << "\t";
                    globals.junk_fileStream << mainMCMC.Qmatrix_gene[i][k];
                    */
                }
                //globals.junk_fileStream << "\n";
            }
            for (int i=0; i<globals.n; i++) {
                for (int k=0; k<K; k++) {
                    Qmatrix_ind_allReps[i][k][mainRep] = mainMCMC.Qmatrix_ind[i][k];
                }
            }
            if (globals.outputQmatrix_pop_on) {
                for (int i=0; i<int(globals.uniquePops.size()); i++) {
                    for (int k=0; k<K; k++) {
                        Qmatrix_pop_allReps[i][k][mainRep] = mainMCMC.Qmatrix_pop[i][k];
                    }
                }
            }
        }
        
        // save harmonic mean of marginal likelihoods
        globals.logEvidence_harmonic[Kindex][mainRep] = mainMCMC.harmonic;
        
        // calculate Structure estimator from joint likelihoods
        globals.structure_loglike_mean[Kindex][mainRep] = mainMCMC.logLikeJoint_sum/double(globals.mainSamples);
        globals.structure_loglike_var[Kindex][mainRep] = mainMCMC.logLikeJoint_sumSquared/double(globals.mainSamples) - globals.structure_loglike_mean[Kindex][mainRep]*globals.structure_loglike_mean[Kindex][mainRep];
        globals.logEvidence_structure[Kindex][mainRep] = globals.structure_loglike_mean[Kindex][mainRep] - 0.5*globals.structure_loglike_var[Kindex][mainRep];
        
        // calculate elements of Evanno's delta K for the *previous* iteration
        if (Kindex>0)
            globals.L_1[Kindex][mainRep] = globals.logEvidence_structure[Kindex][mainRep]-globals.logEvidence_structure[Kindex-1][mainRep];
        if (Kindex>1) {
            globals.L_2[Kindex-1][mainRep] = abs(globals.L_1[Kindex][mainRep]-globals.L_1[Kindex-1][mainRep]);
        }
        
    }
    
    // save Qmatrices
    if (globals.fixLabels_on) {
        for (int i=0; i<globals.geneCopies; i++) {
            for (int k=0; k<K; k++) {
                globals.Qmatrix_gene[Kindex][i][k] = mean(Qmatrix_gene_allReps[i][k]);
                globals.QmatrixError_gene[Kindex][i][k] = sqrt(var(Qmatrix_gene_allReps[i][k])/double(globals.mainRepeats));
            }
        }
        for (int i=0; i<globals.n; i++) {
            for (int k=0; k<K; k++) {
                globals.Qmatrix_ind[Kindex][i][k] = mean(Qmatrix_ind_allReps[i][k]);
                globals.QmatrixError_ind[Kindex][i][k] = sqrt(var(Qmatrix_ind_allReps[i][k]));
            }
        }
        if (globals.outputQmatrix_pop_on) {
            for (int i=0; i<int(globals.uniquePops.size()); i++) {
                for (int k=0; k<K; k++) {
                    globals.Qmatrix_pop[Kindex][i][k] = mean(Qmatrix_pop_allReps[i][k]);
                    globals.QmatrixError_pop[Kindex][i][k] = sqrt(var(Qmatrix_pop_allReps[i][k]));
                }
            }
        }
    }
    
    // calculate grand mean and standard error of harmonic mean estimator
    globals.logEvidence_harmonic_grandMean[Kindex] = mean(globals.logEvidence_harmonic[Kindex]);
    if (globals.mainRepeats>1) {
        globals.logEvidence_harmonic_grandSE[Kindex] = sqrt(var(globals.logEvidence_harmonic[Kindex], 0, true)/globals.mainRepeats);
    }
    
    // calculate grand mean and standard error of Structure estimator
    globals.logEvidence_structure_grandMean[Kindex] = mean(globals.logEvidence_structure[Kindex]);
    if (globals.mainRepeats>1) {
        globals.logEvidence_structure_grandSE[Kindex] = sqrt(var(globals.logEvidence_structure[Kindex], 0, true)/globals.mainRepeats);
    }
    
    // finalise calculation of Evanno's delta K for *previous* iteration
    if (Kindex>1) {
        globals.delta_K[Kindex-1] = mean(globals.L_2[Kindex-1])/sqrt(var(globals.logEvidence_structure[Kindex-1], 0, true));
    }
    
}