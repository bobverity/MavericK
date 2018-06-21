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
MCMC_admixture::MCMC_admixture(globals &globals, int _Kindex) {
    
    // copy some values over from globals object
    uniquePop_counts = globals.uniquePop_counts;
    pop_index = globals.pop_index;
    Kindex = _Kindex;
    K = globals.Kmin+Kindex;
    n = globals.n;
    loci = globals.loci;
    ploidy_vec = globals.ploidy_vec;
    outputQmatrix_pop_on = globals.outputQmatrix_pop_on;
    fixAlpha_on = globals.fixAlpha_on;
    nPops = int(globals.uniquePops.size());
    burnin = globals.burnin;
    samples = globals.samples;
    GTI_pow = globals.GTI_pow;
    
    // create a particle for each rung
    rungs = globals.rungs;
    rung_order = vector<int>(rungs);
    beta_raised_vec = vector<double>(rungs);
    particle_vec = vector<particle_admixture>(rungs);
    for (int rung=0; rung<rungs; rung++) {
        rung_order[rung] = rung;
        beta_raised_vec[rung] = (rungs==1) ? 1 : pow((rung+1)/double(rungs), GTI_pow);
        particle_vec[rung] = particle_admixture(globals, K, globals.alpha[Kindex], beta_raised_vec[rung]);
    }
    coupling_accept = vector<double>(rungs-1);
    
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
    
    // TI point estimates for each rung
    TIpoint_mean = vector<double>(rungs);
    TIpoint_var = vector<double>(rungs);
    TIpoint_SE = vector<double>(rungs);
    
    // overall TI estimate and SE
    logEvidence_TI_store = vector<double>(samples);
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
void MCMC_admixture::perform_MCMC() {
    
    // reset chains
    for (int rung=0; rung<rungs; rung++) {
        particle_vec[rung].reset(true);
    }
    
    // perform MCMC
    for (int rep=0; rep<(burnin+samples); rep++) {
        
        // loop over rungs
        for (int rung=0; rung<rungs; rung++) {
            
            // update group allocation of all individuals
            particle_vec[rung].group_update();
            
            // MH step to update group allocation at individual level. Improves mixing when alpha very small.
            particle_vec[rung].group_update_indLevel();
            
            // MH step to update group allocation at deme level
            //particle_vec[rung].group_update_Klevel();
            
            // if alpha not fixed, update by Metropolis step
            if (fixAlpha_on==0) {
                particle_vec[rung].alpha_update();
            }
            
            // calculate group log-likelihood
            particle_vec[rung].d_logLikeGroup();
        }
        
        // carry out Metropolis coupling
        MetropolisCoupling();
        
        // focus on coldest rung (i.e. the real chain)
        int cold_rung = rung_order[rungs-1];
        
        // fix label-switching problem and save Qmatrix
        updateQmatrix(particle_vec[cold_rung], rep>=burnin);
        
        // draw allele frequencies and calculate joint likelihood
        particle_vec[cold_rung].drawFreqs();
        particle_vec[cold_rung].d_logLikeJoint();
        
        // if no longer in burnin phase
        if (rep>=burnin) {
            
            // add logLikeGroup to running sums
            for (int i=0; i<rungs; i++) {
                logLikeGroup_store[i][rep-burnin] = particle_vec[rung_order[i]].logLikeGroup;
                logLikeGroup_sum[i] += particle_vec[rung_order[i]].logLikeGroup;
                logLikeGroup_sumSquared[i] += particle_vec[rung_order[i]].logLikeGroup*particle_vec[rung_order[i]].logLikeGroup;
            }
            
            // add logLikeJoint to running sums
            logLikeJoint_sum += particle_vec[cold_rung].logLikeJoint;
            logLikeJoint_sumSquared += particle_vec[cold_rung].logLikeJoint*particle_vec[cold_rung].logLikeJoint;
            
            // calculate TI estimate this iteration
            double w1 = pow(1/double(rungs), GTI_pow-1.0);
            logEvidence_TI_store[rep-burnin] = 0.5*GTI_pow*w1*particle_vec[rung_order[0]].logLikeGroup/double(rungs);
            if (rungs>1) {
                for (int i=1; i<rungs; i++) {
                    double loglike1 = particle_vec[rung_order[i-1]].logLikeGroup;
                    double loglike2 = particle_vec[rung_order[i]].logLikeGroup;
                    double w1 = pow(i/double(rungs), GTI_pow-1.0);
                    double w2 = pow((i+1)/double(rungs), GTI_pow-1.0);
                    logEvidence_TI_store[rep-burnin] += 0.5*GTI_pow*(w1*loglike1 + w2*loglike2)/double(rungs);
                }
            }
        }
        
    } // end of MCMC
    
    // finish off acceptance rate vector
    for (int rung=0; rung<int(coupling_accept.size()); rung++) {
        coupling_accept[rung] /= double(burnin+samples);
    }
    //printVector(coupling_accept);
    
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
    
    // calculate TI autocorrelation
    double TI_autoCorr = calculateAutoCorr(logEvidence_TI_store);
    double TI_ESS = samples/TI_autoCorr;
    
    // calculate overall TI point estimate and SE
    logEvidence_TI = mean(logEvidence_TI_store);
    logEvidence_TI_SE = sqrt(var(logEvidence_TI_store)/TI_ESS);
    
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
                Qmatrix_pop[pop_index[i]][k] += Qmatrix_ind[i][k];
            }
        }
        for (int i=0; i<nPops; i++) {
            for (int k=0; k<K; k++) {
                Qmatrix_pop[i][k] /= double(uniquePop_counts[i]);
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
    
    // loop over rungs, starting with the hottest chain and moving to the cold chain. Each time propose a swap with the next rung up.
    for (int i=0; i<(rungs-1); i++) {
        
        // define rungs of interest
        int rung1 = rung_order[i];
        int rung2 = rung_order[i+1];
        
        // get log-likelihoods and beta values of two chains in the comparison
        double loglike1 = particle_vec[rung1].logLikeGroup;
        double loglike2 = particle_vec[rung2].logLikeGroup;
        
        double beta_raised1 = particle_vec[rung1].beta_raised;
        double beta_raised2 = particle_vec[rung2].beta_raised;
        
        // calculate acceptance ratio (still in log space)
        double acceptance = (loglike2*beta_raised1 + loglike1*beta_raised2) - (loglike1*beta_raised1 + loglike2*beta_raised2);
        
        // accept or reject move
        double rand1 = runif1();
        if (log(rand1)<acceptance) {
            
            // swap beta values
            particle_vec[rung1].beta_raised = beta_raised2;
            particle_vec[rung2].beta_raised = beta_raised1;
            
            // swap rung order
            rung_order[i] = rung2;
            rung_order[i+1] = rung1;
            
            // update acceptance rates
            coupling_accept[i]++;
        }
    }
    
}
