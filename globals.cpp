//
//  MavericK
//  globals.cpp
//
//  Created: Bob on 22/09/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "globals.h"

using namespace std;

//------------------------------------------------
// globals::
// constructor for class containing global objects, including parameters and data
globals::globals() {
    
    // file names and paths
	char buffer[255];
    masterRoot_filePath = GETCWD(buffer, sizeof(buffer));
    masterRoot_filePath = masterRoot_filePath + DIRBREAK;
    
    inputRoot_fileName = "";
    outputRoot_fileName = "";
    data_fileName = "data.txt";
    parameters_fileName = "parameters.txt";
    outputLog_fileName = "outputLog.txt";
    outputLikelihood_fileName = "outputLikelihood.csv";
    outputQmatrix_ind_fileName = "outputQmatrix_ind.csv";
    outputQmatrix_pop_fileName = "outputQmatrix_pop.csv";
    outputQmatrix_gene_fileName = "outputQmatrix_gene.csv";
    outputQmatrixError_ind_fileName = "outputQmatrixError_ind.csv";
    outputQmatrixError_pop_fileName = "outputQmatrixError_pop.csv";
    outputQmatrixError_gene_fileName = "outputQmatrixError_gene.csv";
    outputEvidence_fileName = "outputEvidence.csv";
    outputEvidenceNormalised_fileName = "outputEvidenceNormalised.csv";
    outputEvidenceDetails_fileName = "outputEvidenceDetails.csv";
    outputPosteriorGrouping_fileName = "outputPosteriorGrouping.csv";
    outputComparisonStatistics_fileName = "outputComparisonStatistics.csv";
    outputMaxLike_alleleFreqs_fileName = "outputMaxLike_alleleFreqs.csv";
    outputMaxLike_admixFreqs_fileName = "outputMaxLike_admixFreqs.csv";
    
    inputRoot_filePath = masterRoot_filePath + inputRoot_fileName;
    outputRoot_filePath = masterRoot_filePath + outputRoot_fileName;
    data_filePath = inputRoot_filePath + data_fileName;
    parameters_filePath = inputRoot_filePath + parameters_fileName;
    outputLog_filePath = outputRoot_filePath + outputLog_fileName;
    outputLikelihood_filePath = outputRoot_filePath + outputLikelihood_fileName;
    outputQmatrix_ind_filePath = outputRoot_filePath + outputQmatrix_ind_fileName;
    outputQmatrix_pop_filePath = outputRoot_filePath + outputQmatrix_pop_fileName;
    outputQmatrix_gene_filePath = outputRoot_filePath + outputQmatrix_gene_fileName;
    outputQmatrixError_ind_filePath = outputRoot_filePath + outputQmatrixError_ind_fileName;
    outputQmatrixError_pop_filePath = outputRoot_filePath + outputQmatrixError_pop_fileName;
    outputQmatrixError_gene_filePath = outputRoot_filePath + outputQmatrixError_gene_fileName;
    outputEvidence_filePath = outputRoot_filePath + outputEvidence_fileName;
    outputEvidenceNormalised_filePath = outputRoot_filePath + outputEvidenceNormalised_fileName;
    outputEvidenceDetails_filePath = outputRoot_filePath + outputEvidenceDetails_fileName;
    outputComparisonStatistics_filePath = outputRoot_filePath + outputComparisonStatistics_fileName;
    outputMaxLike_alleleFreqs_filePath = outputRoot_filePath + outputMaxLike_alleleFreqs_fileName;
    outputMaxLike_admixFreqs_filePath = outputRoot_filePath + outputMaxLike_admixFreqs_fileName;
    outputPosteriorGrouping_filePath = outputRoot_filePath + outputPosteriorGrouping_fileName;

    // define all default parameter values as pair<string,int> objects, as well as in final class-specific form
    parameterStrings["headerRow_on"] = pair<string,int>("false",0); headerRow_on = false;
    parameterStrings["popCol_on"] = pair<string,int>("false",0); popCol_on = false;
    parameterStrings["ploidyCol_on"] = pair<string,int>("false",0); ploidyCol_on = false;
    parameterStrings["ploidy"] = pair<string,int>("2",0); ploidy = 2;
    parameterStrings["missingData"] = pair<string,int>("-9",0); missingData = "-9";
    
    parameterStrings["Kmin"] = pair<string,int>("1",0); Kmin = 1;
    parameterStrings["Kmax"] = pair<string,int>("2",0); Kmax = 2;
    parameterStrings["admix_on"] = pair<string,int>("false",0); admix_on = false;
    parameterStrings["fixAlpha_on"] = pair<string,int>("true",0); fixAlpha_on = true;
    parameterStrings["alpha"] = pair<string,int>("1.0",0); alpha = vector<double>(1,1.0);
    parameterStrings["alphaPropSD"] = pair<string,int>("0.1",0); vector<double> alphaPropSD(1,0.1);
    
    parameterStrings["exhaustive_on"] = pair<string,int>("false",0); exhaustive_on = false;
    
    parameterStrings["mainBurnin"] = pair<string,int>("100",0); mainBurnin = 100;
    parameterStrings["mainSamples"] = pair<string,int>("1000",0); mainSamples = 1000;
    parameterStrings["mainRungs"] = pair<string,int>("10",0); mainRungs = 10;
    
    parameterStrings["EMalgorithm_on"] = pair<string,int>("false",0); EMalgorithm_on = false;
    parameterStrings["EMrepeats"] = pair<string,int>("10",0); EMrepeats = 10;
    parameterStrings["EMiterations"] = pair<string,int>("100",0); EMiterations = 100;
    
    parameterStrings["outputLog_on"] = pair<string,int>("true",0); outputLog_on = true;
    parameterStrings["outputLikelihood_on"] = pair<string,int>("false",0); outputLikelihood_on = false;
    parameterStrings["outputQmatrix_ind_on"] = pair<string,int>("true",0); outputQmatrix_ind_on = true;
    parameterStrings["outputQmatrix_pop_on"] = pair<string,int>("false",0); outputQmatrix_pop_on = false;
    parameterStrings["outputQmatrix_gene_on"] = pair<string,int>("false",0); outputQmatrixError_gene_on = false;
    parameterStrings["outputQmatrixError_ind_on"] = pair<string,int>("false",0); outputQmatrixError_ind_on = false;
    parameterStrings["outputQmatrixError_pop_on"] = pair<string,int>("false",0); outputQmatrixError_pop_on = false;
    parameterStrings["outputQmatrixError_gene_on"] = pair<string,int>("false",0); outputQmatrixError_gene_on = false;
    parameterStrings["outputEvidence_on"] = pair<string,int>("true",0); outputEvidence_on = true;
    parameterStrings["outputEvidenceNormalised_on"] = pair<string,int>("true",0); outputEvidenceNormalised_on = true;
    parameterStrings["outputEvidenceDetails_on"] = pair<string,int>("false",0); outputEvidenceDetails_on = false;
    parameterStrings["outputPosteriorGrouping_on"] = pair<string,int>("false",0); outputPosteriorGrouping_on = false;
    parameterStrings["outputComparisonStatistics_on"] = pair<string,int>("false",0); outputComparisonStatistics_on = false;
    parameterStrings["outputMaxLike_alleleFreqs_on"] = pair<string,int>("false",0); outputMaxLike_alleleFreqs_on = false;
    parameterStrings["outputMaxLike_admixFreqs_on"] = pair<string,int>("false",0); outputMaxLike_admixFreqs_on = false;
    
    parameterStrings["outputQmatrix_structureFormat_on"] = pair<string,int>("false",0); outputQmatrix_structureFormat_on = false;
    parameterStrings["suppressWarning1_on"] = pair<string,int>("false",0); suppressWarning1_on = false;
    
    // parameters not defined by user
    lambda = 1.0;
};

//------------------------------------------------
// initialise global objects with empty values
void initialiseGlobals(globals &globals) {
    
    // Qmatrices
    globals.Qmatrix_gene = vector< vector< vector<double> > >(globals.Kmax-globals.Kmin+1);
    globals.QmatrixError_gene = vector< vector< vector<double> > >(globals.Kmax-globals.Kmin+1);
    globals.Qmatrix_ind = vector< vector< vector<double> > >(globals.Kmax-globals.Kmin+1);
    globals.QmatrixError_ind = vector< vector< vector<double> > >(globals.Kmax-globals.Kmin+1);
    globals.Qmatrix_pop = vector< vector< vector<double> > >(globals.Kmax-globals.Kmin+1);
    globals.QmatrixError_pop = vector< vector< vector<double> > >(globals.Kmax-globals.Kmin+1);
    
    for (int Kindex=0; Kindex<(globals.Kmax-globals.Kmin+1); Kindex++) {
        int K = globals.Kmin+Kindex;
        globals.Qmatrix_gene[Kindex] = vector< vector<double> >(globals.geneCopies,vector<double>(K));
        globals.QmatrixError_gene[Kindex] = vector< vector<double> >(globals.geneCopies,vector<double>(K));
        globals.Qmatrix_ind[Kindex] = vector< vector<double> >(globals.n,vector<double>(K));
        globals.QmatrixError_ind[Kindex] = vector< vector<double> >(globals.n,vector<double>(K));
        globals.Qmatrix_pop[Kindex] = vector< vector<double> >(globals.uniquePops.size(),vector<double>(K));
        globals.QmatrixError_pop[Kindex] = vector< vector<double> >(globals.uniquePops.size(),vector<double>(K));
    }
    
    // evidence estimates
    vector<double> nanVec = vector<double>(globals.Kmax-globals.Kmin+1,-sqrt(-1.0));
    globals.logEvidence_exhaustive = nanVec;
    globals.logEvidence_harmonic = vector<double>(globals.Kmax-globals.Kmin+1);
    globals.structure_loglike_mean = vector<double>(globals.Kmax-globals.Kmin+1);
    globals.structure_loglike_var = vector<double>(globals.Kmax-globals.Kmin+1);
    globals.logEvidence_structure = vector<double>(globals.Kmax-globals.Kmin+1);
    
    globals.TIpoint_mean = vector< vector<double> >(globals.Kmax-globals.Kmin+1,vector<double>(globals.mainRungs,-sqrt(-1.0)));
    globals.TIpoint_var = vector< vector<double> >(globals.Kmax-globals.Kmin+1,vector<double>(globals.mainRungs,-sqrt(-1.0)));
    globals.TIpoint_SE = vector< vector<double> >(globals.Kmax-globals.Kmin+1,vector<double>(globals.mainRungs,-sqrt(-1.0)));
    globals.logEvidence_TI = nanVec;
    globals.logEvidence_TI_SE = nanVec;
    
    globals.posterior_exhaustive = nanVec;
    globals.posterior_harmonic_mean = nanVec;
    globals.posterior_harmonic_LL = nanVec;
    globals.posterior_harmonic_UL = nanVec;
    globals.posterior_structure_mean = nanVec;
    globals.posterior_structure_LL = nanVec;
    globals.posterior_structure_UL = nanVec;
    globals.posterior_TI_mean = nanVec;
    globals.posterior_TI_LL = nanVec;
    globals.posterior_TI_UL = nanVec;
    
    globals.AIC = nanVec;
    globals.BIC = nanVec;
    globals.DIC_Gelman = nanVec;
    globals.DIC_Spiegelhalter = nanVec;
    
}
