//
//  MavericK
//  TI.cpp
//
//  Created: Bob on 23/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "TI.h"

using namespace std;

//------------------------------------------------
// thermodynamic integral estimator for no-admixture model
void TI_noAdmixture(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // special case if K==1
    if (K==1) {
        for (int TIrep=0; TIrep<globals.thermodynamicRungs; TIrep++) {
            globals.TIpoint_mean[Kindex][TIrep] = globals.logEvidence_exhaustive[Kindex];
            globals.TIpoint_var[Kindex][TIrep] = 0;
            globals.TIpoint_SE[Kindex][TIrep] = 0;
        }
        globals.logEvidence_TI[Kindex] = globals.logEvidence_exhaustive[Kindex];
        globals.logEvidence_TI_SE[Kindex] = 0;
        return;
    }
    
    // set up beta vector
    vector<double> betaVec(globals.thermodynamicRungs);
    for (int TIrep=0; TIrep<globals.thermodynamicRungs; TIrep++)
        betaVec[TIrep] = double(TIrep)/(globals.thermodynamicRungs-1);
    
    // carry out thermodynamic integration
    for (int TIrep=0; TIrep<globals.thermodynamicRungs; TIrep++) {
        double beta = betaVec[TIrep];
        
        char * buffer = new char[255];
        sprintf(buffer, "%.2f", beta);
        string s = buffer;
        coutAndLog("  power = "+s+"\n", globals.outputLog_on, globals.outputLog_fileStream);
        
        // define MCMC object
        MCMCobject_noAdmixture TI_MCMC(globals, Kindex, globals.thermodynamicBurnin, globals.thermodynamicSamples, globals.thermodynamicThinning, beta);
        
        // perform MCMC
        TI_MCMC.reset(true);
        TI_MCMC.perform_MCMC(globals, false, true, false, false, false, 1);
        
        // calculate autocorrelation
        double autoCorr = calculateAutoCorr(TI_MCMC.logLikeGroup_store);
        double ESS = globals.thermodynamicSamples/autoCorr;
        
        // save results of this power
        globals.TIpoint_mean[Kindex][TIrep] = TI_MCMC.logLikeGroup_sum/double(globals.thermodynamicSamples);
        globals.TIpoint_var[Kindex][TIrep] = TI_MCMC.logLikeGroup_sumSquared/double(globals.thermodynamicSamples) - globals.TIpoint_mean[Kindex][TIrep]*globals.TIpoint_mean[Kindex][TIrep];
		// avoid arithmetic underflow
		if (globals.TIpoint_var[Kindex][TIrep]<0)
			globals.TIpoint_var[Kindex][TIrep] = 0;
        globals.TIpoint_SE[Kindex][TIrep] = sqrt(globals.TIpoint_var[Kindex][TIrep]/ESS);
        
    }
    
    // calculate thermodynamic integral estimate. Note that at this stage we assume equally spaced rungs.
    double Dsum = 0.5*globals.TIpoint_mean[Kindex][0] + 0.5*globals.TIpoint_mean[Kindex][globals.thermodynamicRungs-1];
    double Vsum = 0.25*globals.TIpoint_SE[Kindex][0]*globals.TIpoint_SE[Kindex][0] + 0.25*globals.TIpoint_SE[Kindex][globals.thermodynamicRungs-1]*globals.TIpoint_SE[Kindex][globals.thermodynamicRungs-1];
    if (globals.thermodynamicRungs>2) {
        for (int TIrep=1; TIrep<(globals.thermodynamicRungs-1); TIrep++) {
            Dsum += globals.TIpoint_mean[Kindex][TIrep];
            Vsum += globals.TIpoint_SE[Kindex][TIrep]*globals.TIpoint_SE[Kindex][TIrep];
        }
    }
    Dsum /= double(globals.thermodynamicRungs-1);
    Vsum /= double(globals.thermodynamicRungs-1)*double(globals.thermodynamicRungs-1);
    globals.logEvidence_TI[Kindex] = Dsum;
    globals.logEvidence_TI_SE[Kindex] = sqrt(Vsum);
    
}

//------------------------------------------------
// thermodynamic integral estimator for admixture model
void TI_admixture(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // special case if K==1
    if (K==1) {
        for (int TIrep=0; TIrep<globals.thermodynamicRungs; TIrep++) {
            globals.TIpoint_mean[Kindex][TIrep] = globals.logEvidence_exhaustive[Kindex];
            globals.TIpoint_var[Kindex][TIrep] = 0;
            globals.TIpoint_SE[Kindex][TIrep] = 0;
        }
        globals.logEvidence_TI[Kindex] = globals.logEvidence_exhaustive[Kindex];
        globals.logEvidence_TI_SE[Kindex] = 0;
        return;
    }
    
    // set up beta vector
    vector<double> betaVec(globals.thermodynamicRungs);
    for (int TIrep=0; TIrep<globals.thermodynamicRungs; TIrep++)
        betaVec[TIrep] = double(TIrep)/(globals.thermodynamicRungs-1);
    
    // carry out thermodynamic integration
    for (int TIrep=0; TIrep<globals.thermodynamicRungs; TIrep++) {
        double beta = betaVec[TIrep];
        
        char * buffer = new char[8];
        sprintf(buffer, "%.2f", beta);
        string s = buffer;
        coutAndLog("  power = "+s+"\n", globals.outputLog_on, globals.outputLog_fileStream);
        
        // define MCMC object
        MCMCobject_admixture TI_MCMC(globals, Kindex, globals.thermodynamicBurnin, globals.thermodynamicSamples, globals.thermodynamicThinning, beta);
        
        // perform MCMC
        TI_MCMC.reset(true);
		TI_MCMC.perform_MCMC(globals, false, true, false, false, false, TIrep);
        
        // calculate autocorrelation
        double autoCorr = calculateAutoCorr(TI_MCMC.logLikeGroup_store);
        double ESS = globals.thermodynamicSamples/autoCorr;
        
        // save results of this power
        globals.TIpoint_mean[Kindex][TIrep] = TI_MCMC.logLikeGroup_sum/double(globals.thermodynamicSamples);
        globals.TIpoint_var[Kindex][TIrep] = TI_MCMC.logLikeGroup_sumSquared/double(globals.thermodynamicSamples) - globals.TIpoint_mean[Kindex][TIrep]*globals.TIpoint_mean[Kindex][TIrep];
        globals.TIpoint_SE[Kindex][TIrep] = sqrt(globals.TIpoint_var[Kindex][TIrep]/ESS);
        
    }
    
    // calculate thermodynamic integral estimate. Note that at this stage we assume equally spaced rungs.
    double Dsum = 0.5*globals.TIpoint_mean[Kindex][0] + 0.5*globals.TIpoint_mean[Kindex][globals.thermodynamicRungs-1];
    double Vsum = 0.25*globals.TIpoint_SE[Kindex][0]*globals.TIpoint_SE[Kindex][0] + 0.25*globals.TIpoint_SE[Kindex][globals.thermodynamicRungs-1]*globals.TIpoint_SE[Kindex][globals.thermodynamicRungs-1];
    if (globals.thermodynamicRungs>2) {
        for (int TIrep=1; TIrep<(globals.thermodynamicRungs-1); TIrep++) {
            Dsum += globals.TIpoint_mean[Kindex][TIrep];
            Vsum += globals.TIpoint_SE[Kindex][TIrep]*globals.TIpoint_SE[Kindex][TIrep];
        }
    }
    Dsum /= double(globals.thermodynamicRungs-1);
    Vsum /= double(globals.thermodynamicRungs-1)*double(globals.thermodynamicRungs-1);
    globals.logEvidence_TI[Kindex] = Dsum;
    globals.logEvidence_TI_SE[Kindex] = sqrt(Vsum);
    
}