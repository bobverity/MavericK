//
//  MavericK
//  MCMC_admixture.cpp
//
//  Created: Bob on 06/11/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "MCMC_admixture.h"

using namespace std;

//------------------------------------------------
// MCMC_admixture::
// constructor for class containing all elements required for MCMC under admixture model
MCMC_admixture::MCMC_admixture(globals &globals, int _Kindex, int _rungs) {
    
    // copy some values over from globals object
    Kindex = _Kindex;
    K = globals.Kmin+Kindex;
    n = globals.n;
    loci = globals.loci;
    ploidy_vec = globals.ploidy_vec;
    outputQmatrix_pop_on = globals.outputQmatrix_pop_on;
    nPops = int(globals.uniquePops.size());
    burnin = globals.mainBurnin;
    samples = globals.mainSamples;
    
    // create a particle for each rung
    rungs = _rungs;
    rungOrder = vector<int>(rungs);
    betaVec = vector<double>(rungs);
    particleVec = vector<particle_admixture>(rungs);
    for (int rung=0; rung<rungs; rung++) {
        rungOrder[rung] = rung;
        betaVec[rung] = rung/double(rungs-1);
        particleVec[rung] = particle_admixture(globals, K, globals.alpha[Kindex], globals.alphaPropSD[Kindex], betaVec[rung]);
    }
    acceptanceRate = vector<double>(rungs-1);
    
    // likelihoods for each rung
    logLikeGroup_sum = vector<double>(rungs);
    logLikeGroup_sumSquared = vector<double>(rungs);
    logLikeGroup_store = vector< vector<double> >(rungs, vector<double>(samples));
    logLikeJoint_sum = 0;
    logLikeJoint_sumSquared = 0;
    
    // initialise Qmatrices
    logQmatrix_gene_running = vector< vector< vector< vector<double> > > >(n, vector< vector< vector<double> > >(loci));
    for (int ind=0; ind<n; ind++) {
        logQmatrix_gene_running[ind] = vector< vector< vector<double> > >(loci, vector< vector<double> >(ploidy_vec[ind], vector<double>(K)));
    }
    logQmatrix_gene_update = logQmatrix_gene_running;
    Qmatrix_gene_update = logQmatrix_gene_running;
    
    Qmatrix_gene = logQmatrix_gene_running;
    Qmatrix_ind = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_pop = vector< vector<double> >(nPops, vector<double>(K));
    
    // initialise ordering of labels
    labelOrder = vector<int>(K);
    labelOrder_new = vector<int>(K);
    for (int k=0; k<K; k++) {
        labelOrder[k] = k;
    }
    
    // autocorrelation for each rung
    autoCorr = vector<double>(rungs);
    ESS = vector<double>(rungs);
    
    // harmonic mean estimator
    harmonic = NEGINF;
    
    // TI point estimates for each rung
    TIpoint_mean = vector<double>(rungs);
    TIpoint_var = vector<double>(rungs);
    TIpoint_SE = vector<double>(rungs);
    
    // overall TI point estimate and SE
    logEvidence_TI = 0;
    logEvidence_TI_var = 0;
    logEvidence_TI_SE = 0;
    
    // initialise objects for Hungarian algorithm
    costMat = vector< vector<double> >(K, vector<double>(K));
    bestPermOrder = vector<int>(K);
    
    edgesLeft = vector<int>(K);
    edgesRight = vector<int>(K);
    blockedLeft = vector<int>(K);
    blockedRight = vector<int>(K);
    
}

//------------------------------------------------
// MCMC_admixture::
// perform complete MCMC under admixture model
void MCMC_admixture::perform_MCMC(globals &globals) {
    
    // reset chains
    for (int rung=0; rung<rungs; rung++) {
        particleVec[rung].reset(true);
    }
    
    // perform MCMC
    for (int rep=0; rep<(burnin+samples); rep++) {
        
        // loop over rungs
        for (int rung=0; rung<rungs; rung++) {
            
            // update group allocation of all individuals
            particleVec[rung].group_update();
            
            // MH step to update group allocation at individual level. Improves mixing when alpha very small.
            particleVec[rung].group_update_indLevel();
            
            // MH step to update group allocation at deme level
            //particleVec[rung].group_update_Klevel();
            
            // if alpha not fixed, update by Metropolis step
            if (globals.fixAlpha_on==0)
               particleVec[rung].alpha_update();
            
            // calculate group log-likelihood
            particleVec[rung].d_logLikeGroup();
        }
        
        // carry out Metropolis coupling
        MetropolisCoupling();
        
        // PRINT JUNK TO FILE
        for (int i=0; i<rungs; i++) {
            globals.junk_fileStream << particleVec[rungOrder[i]].logLikeGroup << "\t";
            //globals.junk_fileStream << particleVec[rungOrder[i]].alpha << "\t";
        }
        globals.junk_fileStream << "\n";
        
        // focus on coldest rung (i.e. the real chain)
        rung1 = rungOrder[rungs-1];
        
        //print(particleVec[rung1].alpha);
        
        // fix label-switching problem and save Qmatrix
        updateQmatrix(particleVec[rung1], rep>=burnin);
        
        // if no longer in burnin phase
        if (rep>=burnin) {
            
            // add logLikeGroup to running sums
            for (int i=0; i<rungs; i++) {
                logLikeGroup_store[i][rep-burnin] = particleVec[rungOrder[i]].logLikeGroup;
                logLikeGroup_sum[i] += particleVec[rungOrder[i]].logLikeGroup;
                logLikeGroup_sumSquared[i] += particleVec[rungOrder[i]].logLikeGroup*particleVec[rungOrder[i]].logLikeGroup;
            }
            
            // add logLikeGroup to harmonic mean estimator
            harmonic = logSum(harmonic, -particleVec[rung1].logLikeGroup);
            
            // draw allele frequencies and calculate joint likelihood
            particleVec[rung1].drawFreqs();
            particleVec[rung1].d_logLikeJoint();
            logLikeJoint_sum += particleVec[rung1].logLikeJoint;
            logLikeJoint_sumSquared += particleVec[rung1].logLikeJoint*particleVec[rung1].logLikeJoint;
            
        }
        
    } // end of MCMC
    
    // finish off acceptance rate vector
    for (int rung=0; rung<rungs; rung++) {
        acceptanceRate[rung] /= double(burnin+samples);
    }
    printVector(acceptanceRate);
    
    // finish off harmonic mean estimator
    harmonic = log(double(samples))-harmonic;
    
    // calculate Structure estimator from joint likelihoods
    structure_loglike_mean = logLikeJoint_sum/double(samples);
    structure_loglike_var = logLikeJoint_sumSquared/double(samples) - structure_loglike_mean*structure_loglike_mean;
    logEvidence_structure = structure_loglike_mean - 0.5*structure_loglike_var;
    
    // process logLikeGroup values for each rung. Calculate mean and SE, taking into account autocorrelation.
    for (int rung=0; rung<rungs; rung++) {
        
        // calculate autocorrelation
        autoCorr[rung] = calculateAutoCorr(logLikeGroup_store[rung]);
        ESS[rung] = samples/autoCorr[rung];
        
        // calculate mean and SE of log-likelihood
        TIpoint_mean[rung] = logLikeGroup_sum[rung]/double(samples);
        TIpoint_var[rung] = logLikeGroup_sumSquared[rung]/double(samples) - TIpoint_mean[rung]*TIpoint_mean[rung];
        if (TIpoint_var[rung]<0) {   // avoid arithmetic underflow
            TIpoint_var[rung] = 0;
        }
        TIpoint_SE[rung] = sqrt(TIpoint_var[rung]/ESS[rung]);
    }
    
    // calculate overall TI point estimate and SE. Note that simply summing the variance by looping over rungs will give the wrong answer, as it fails to account for the fact that variance(2*A) is not 2*variance(A) but rather 4*variance(A).
    logEvidence_TI = 0.5*TIpoint_mean[0]*(betaVec[1]-betaVec[0]) + 0.5*TIpoint_mean[rungs-1]*(betaVec[rungs-1]-betaVec[rungs-2]);
    logEvidence_TI_var = 0.25*TIpoint_var[0]/ESS[0]*pow(betaVec[1]-betaVec[0],2) + 0.25*TIpoint_var[rungs-1]/ESS[rungs-1]*pow(betaVec[rungs-1]-betaVec[rungs-2],2);
    if (rungs>2) {
        for (int rung=1; rung<(rungs-1); rung++) {
            logEvidence_TI += TIpoint_mean[rung]*(betaVec[rung] - betaVec[rung-1]);
            logEvidence_TI_var += TIpoint_var[rung]/ESS[rung]*pow(betaVec[rung]-betaVec[rung-1],2);
        }
    }
    logEvidence_TI_SE = sqrt(logEvidence_TI_var);
    
    // calculate individual level Qmatrices
    for (int ind=0; ind<n; ind++) {
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                for (int k=0; k<K; k++) {
                    Qmatrix_ind[ind][k] += Qmatrix_gene[ind][l][p][k];
                }
            }
        }
        for (int k=0; k<K; k++) {
            Qmatrix_ind[ind][k] /= double(ploidy_vec[ind]*loci);
        }
    }
    
    // calculate population level Qmatrices
    if (outputQmatrix_pop_on) {
        for (int i=0; i<n; i++) {
            for (int k=0; k<K; k++) {
                Qmatrix_pop[globals.pop_index[i]][k] += Qmatrix_ind[i][k];
            }
        }
        for (int i=0; i<nPops; i++) {
            for (int k=0; k<K; k++) {
                Qmatrix_pop[i][k] /= double(globals.uniquePop_counts[i]);
            }
        }
    }
    
}

//------------------------------------------------
// MCMC_admixture::
// choose best permutation of labels using method of Stephens (2000)
void MCMC_admixture::updateQmatrix(particle_admixture &particle, bool outOfBurnin) {
    
    // calculate update to Qmatrix
    for (int ind=0; ind<n; ind++) {
        for (int l=0; l<loci; l++) {
            for (int p=0; p<particle.ploidy_vec[ind]; p++) {
                particle.group_probs(ind,l,p);
                for (int k=0; k<K; k++) {
                    logQmatrix_gene_update[ind][l][p][k] = particle.logProbVec[k];
                    Qmatrix_gene_update[ind][l][p][k] = particle.probVec[k]/particle.probVecSum;
                }
            }
        }
    }
    
    // calculate cost matrix from old and new Qmatrices
    for (int k1=0; k1<K; k1++) {
        fill(costMat[k1].begin(), costMat[k1].end(), 0);
        for (int k2=0; k2<K; k2++) {
            for (int ind=0; ind<n; ind++) {
                for (int l=0; l<loci; l++) {
                    for (int p=0; p<ploidy_vec[ind]; p++) {
                        costMat[k1][k2] += Qmatrix_gene_update[ind][l][p][labelOrder[k1]]*(logQmatrix_gene_update[ind][l][p][labelOrder[k1]]-logQmatrix_gene_running[ind][l][p][k2]);
                    }
                }
            }
        }
    }
    
    // find best permutation of current labels
    bestPerm = hungarian(costMat, edgesLeft, edgesRight, blockedLeft, blockedRight);
    
    // define bestPermOrder. If the numbers 1:m_K are placed in best-perm-order then we arrive back at bestPerm. In R terms we would say bestPermOrder=order(bestPerm).
    for (int k=0; k<K; k++) {
        bestPermOrder[bestPerm[k]] = k;
    }
    
    // replace old label order with new
    for (int k=0; k<K; k++) {
        labelOrder_new[k] = labelOrder[bestPermOrder[k]];
    }
    labelOrder = labelOrder_new;
    
    // add logQmatrix_ind_new to logQmatrix_gene_running
    for (int ind=0; ind<n; ind++) {
        for (int l=0; l<loci; l++) {
            for (int p=0; p<ploidy_vec[ind]; p++) {
                for (int k=0; k<K; k++) {
                    logQmatrix_gene_running[ind][l][p][k] = logSum(logQmatrix_gene_running[ind][l][p][k], logQmatrix_gene_update[ind][l][p][labelOrder[k]]);
                }
            }
        }
    }
    
    // store Qmatrix values
    if (outOfBurnin) {
        for (int ind=0; ind<n; ind++) {
            for (int l=0; l<loci; l++) {
                for (int p=0; p<ploidy_vec[ind]; p++) {
                    for (int k=0; k<K; k++) {
                        Qmatrix_gene[ind][l][p][k] += Qmatrix_gene_update[ind][l][p][labelOrder[k]]/double(samples);
                    }
                }
            }
        }
    }
    
}

//------------------------------------------------
// MCMC_admixture::
// swap chains in Metropolis step
void MCMC_admixture::MetropolisCoupling() {
    
    // loop over rungs, starting with the hottest chain and moving to the cold chain. Each time propose a swap with a randomly chosen chain.
    for (int i1=0; i1<rungs; i1++) {
        
        // draw value i2!=i1
        int i2 = sample2(1,rungs)-1;
        if (i2==i1) {
            i2++;
            if (i2==rungs)
                i2 = 0;
        }
        rung1 = rungOrder[i1];
        rung2 = rungOrder[i2];
        
        // get log-likelihoods of two chains in the comparison
        double logLike1 = particleVec[rung1].logLikeGroup;
        double logLike2 = particleVec[rung2].logLikeGroup;
        
        // calculate acceptance ratio (still in log space)
        double acceptance = (logLike2*betaVec[rung1] + logLike1*betaVec[rung2]) - (logLike1*betaVec[rung1] + logLike2*betaVec[rung2]);
        
        // accept or reject move
        double rand1 = runif1();
        if (log(rand1)<acceptance) {
            
            // swap beta values
            double spareBeta = particleVec[rung1].beta;
            particleVec[rung1].beta = particleVec[rung2].beta;
            particleVec[rung2].beta = spareBeta;
            
            // swap rung order
            rungOrder[i1] = rung2;
            rungOrder[i2] = rung1;
            
            // update Metropolis coupling acceptance ratio
            acceptanceRate[i1] ++;
            
        }
    }
    
}
