//
//  MavericK
//  globals.h
//
//  Created: Bob on 22/09/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Defines a class containing data and parameters that will be used throughout the program.
//
// ---------------------------------------------------------------------------

#include <random>
#include <fstream>
#include <sstream>
#include <map>
#include "OSfunctions.h"

#ifndef __Maverick1_0__globals__
#define __Maverick1_0__globals__

//------------------------------------------------
// class containing global objects, including parameters and data
class globals {
    
public:
    
    // PUBLIC OBJECTS
    
    // file names and paths
    std::string masterRoot_filePath;
    
    std::string inputRoot_fileName;
    std::string outputRoot_fileName;
    std::string data_fileName;
    std::string parameters_fileName;
    std::string outputLog_fileName;
    std::string outputLikelihood_fileName;
    std::string outputQmatrix_ind_fileName;
    std::string outputQmatrix_pop_fileName;
    std::string outputQmatrix_gene_fileName;
    std::string outputQmatrixError_ind_fileName;
    std::string outputQmatrixError_pop_fileName;
    std::string outputQmatrixError_gene_fileName;
    std::string outputEvidence_fileName;
    std::string outputEvidenceNormalised_fileName;
    std::string outputEvidenceDetails_fileName;
    std::string outputPosteriorGrouping_fileName;
    std::string outputComparisonStatistics_fileName;
    std::string outputEvanno_fileName;
    std::string outputMaxLike_alleleFreqs_fileName;
    std::string outputMaxLike_admixFreqs_fileName;
    
    std::string inputRoot_filePath;
    std::string outputRoot_filePath;
    std::string data_filePath;
    std::string parameters_filePath;
    std::string outputLog_filePath;
    std::string outputLikelihood_filePath;
    std::string outputQmatrix_ind_filePath;
    std::string outputQmatrix_pop_filePath;
    std::string outputQmatrix_gene_filePath;
    std::string outputQmatrixError_ind_filePath;
    std::string outputQmatrixError_pop_filePath;
    std::string outputQmatrixError_gene_filePath;
    std::string outputAdmixture_filePath;
    std::string outputEvidence_filePath;
    std::string outputEvidenceNormalised_filePath;
    std::string outputEvidenceDetails_filePath;
    std::string outputPosteriorGrouping_filePath;
    std::string outputComparisonStatistics_filePath;
    std::string outputEvanno_filePath;
    std::string outputMaxLike_alleleFreqs_filePath;
    std::string outputMaxLike_admixFreqs_filePath;
    std::string junk_filePath;
    
    // file streams
    std::ifstream parameters_fileStream;
    std::ifstream data_fileStream;
    std::ofstream outputLog_fileStream;
    std::ofstream outputLikelihood_fileStream;
    std::ofstream outputQmatrix_ind_fileStream;
    std::ofstream outputQmatrix_pop_fileStream;
    std::ofstream outputQmatrix_gene_fileStream;
    std::ofstream outputQmatrixError_ind_fileStream;
    std::ofstream outputQmatrixError_pop_fileStream;
    std::ofstream outputQmatrixError_gene_fileStream;
    std::ofstream outputEvidence_fileStream;
    std::ofstream outputEvidenceNormalised_fileStream;
    std::ofstream outputEvidenceDetails_fileStream;
    std::ofstream outputPosteriorGrouping_fileStream;
    std::ofstream outputComparisonStatistics_fileStream;
    std::ofstream outputEvanno_fileStream;
    std::ofstream outputMaxLike_alleleFreqs_fileStream;
    std::ofstream outputMaxLike_admixFreqs_fileStream;
    std::ofstream junk_fileStream;
    
    
    // parameters from file
    std::map< std::string, std::pair<std::string,int> > parameterStrings;
    
    bool headerRow_on;
    bool popCol_on;
    bool ploidyCol_on;
    int ploidy;
    std::string missingData;
    
    int Kmin;
    int Kmax;
    bool admix_on;
    bool fixAlpha_on;
    std::vector<double> alpha;
    std::vector<double> alphaPropSD;
    
    bool exhaustive_on;
    int mainRepeats;
    int mainBurnin;
    int mainSamples;
    int mainThinning;
    bool thermodynamic_on;
    int thermodynamicRungs;
    int thermodynamicBurnin;
    int thermodynamicSamples;
    int thermodynamicThinning;
    bool EMalgorithm_on;
    int EMrepeats;
    int EMiterations;
    
    bool outputLog_on;
    bool outputLikelihood_on;
    bool outputQmatrix_ind_on;
    bool outputQmatrix_pop_on;
    bool outputQmatrix_gene_on;
    bool outputQmatrixError_ind_on;
    bool outputQmatrixError_pop_on;
    bool outputQmatrixError_gene_on;
    bool outputAdmixture_on;
    bool outputEvidence_on;
    bool outputEvidenceNormalised_on;
    bool outputEvidenceDetails_on;
    bool outputPosteriorGrouping_on;
    bool outputComparisonStatistics_on;
    bool outputEvanno_on;
    bool outputMaxLike_alleleFreqs_on;
    bool outputMaxLike_admixFreqs_on;
    
    bool outputQmatrix_structureFormat_on;
    bool suppressWarning1_on;
    bool fixLabels_on;
    
    // parameters not defined by user
    double lambda;
    
    
    // data
    std::vector<std::string> indLabels_vec;
    std::vector<std::string> pop_vec;
    std::vector<std::string> uniquePops;
    std::vector<int> uniquePop_counts;
    std::vector<int> pop_index;
    std::vector<int> ploidy_vec;
    std::vector<int> missing_vec;
    std::vector< std::vector< std::vector<int> > > data;
    int n;
    int loci;
    std::vector<int> J;
    std::vector< std::vector<std::string> > uniqueAlleles;
    int geneCopies;
    
    // log lookup table
    std::vector< std::vector<double> > log_lookup;
    
    // objects for storing results
    std::vector< std::vector< std::vector<double> > > Qmatrix_gene;
    std::vector< std::vector< std::vector<double> > > QmatrixError_gene;
    std::vector< std::vector< std::vector<double> > > Qmatrix_ind;
    std::vector< std::vector< std::vector<double> > > QmatrixError_ind;
    std::vector< std::vector< std::vector<double> > > Qmatrix_pop;
    std::vector< std::vector< std::vector<double> > > QmatrixError_pop;
    
    std::vector<double> logEvidence_exhaustive;
    
    std::vector< std::vector<double> > logEvidence_harmonic;
    std::vector<double> logEvidence_harmonic_grandMean;
    std::vector<double> logEvidence_harmonic_grandSE;
    
    std::vector< std::vector<double> > structure_loglike_mean;
    std::vector< std::vector<double> > structure_loglike_var;
    std::vector< std::vector<double> > logEvidence_structure;
    std::vector<double> logEvidence_structure_grandMean;
    std::vector<double> logEvidence_structure_grandSE;
    
    std::vector< std::vector<double> > TIpoint_mean;
    std::vector< std::vector<double> > TIpoint_var;
    std::vector< std::vector<double> > TIpoint_SE;
    std::vector<double> logEvidence_TI;
    std::vector<double> logEvidence_TI_SE;
    
    std::vector<double> posterior_exhaustive;
    std::vector<double> posterior_harmonic_mean;
    std::vector<double> posterior_harmonic_LL;
    std::vector<double> posterior_harmonic_UL;
    std::vector<double> posterior_structure_mean;
    std::vector<double> posterior_structure_LL;
    std::vector<double> posterior_structure_UL;
    std::vector<double> posterior_TI_mean;
    std::vector<double> posterior_TI_LL;
    std::vector<double> posterior_TI_UL;
    
    double maxLike;
    std::vector< std::vector< std::vector<double> > > max_alleleFreqs;
    std::vector< std::vector<double> > max_admixFreqs;
    
    std::vector<double> AIC;
    std::vector<double> BIC;
    std::vector<double> DIC_Spiegelhalter;
    std::vector<double> DIC_Gelman;
    
    // first and second derivative of L (structure estimator) for Evanno calculation
    std::vector< std::vector<double> > L_1;
    std::vector< std::vector<double> > L_2;
    std::vector<double> delta_K;
    
    // PUBLIC FUNCTIONS
    
    // core functions
    globals(); // constructor
    
};

#endif
