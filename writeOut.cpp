//
//  MavericK
//  writeOut.cpp
//
//  Created: Bob on 25/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "writeOut.h"

using namespace std;

//------------------------------------------------
// safely open output file stream, otherwise error. Option to write errors to file if writeErrorToFile is true.
ofstream safe_ofstream(string fileName, bool writeToFile, ofstream &logFileStream) {
    ofstream fileStream(fileName);
    if (!fileStream.is_open()) {
        cerrAndLog("\nError: failed to write to file: "+fileName+string("\n"), writeToFile, logFileStream);
        exit(1);
    }
    return(fileStream);
}

//------------------------------------------------
// open file streams that are common to all K
void openFileStreams(globals &globals) {
    
    // Likelihood
    if (globals.outputLikelihood_on) {
        // open file stream
        globals.outputLikelihood_fileStream = safe_ofstream(globals.outputLikelihood_filePath, globals.outputLog_on, globals.outputLog_fileStream);
        // fixed headers
        globals.outputLikelihood_fileStream << "K,mainRep,MCMCsample,loglike_marginal,loglike_joint";
        // alpha if admixture model
        if (globals.admix_on)
            globals.outputLikelihood_fileStream << ",alpha";
        
        globals.outputLikelihood_fileStream << "\n";
        globals.outputLikelihood_fileStream.flush();
    }
    
    // Evidence
    if (globals.outputEvidence_on) {
        // open file stream
        globals.outputEvidence_fileStream = safe_ofstream(globals.outputEvidence_filePath, globals.outputLog_on, globals.outputLog_fileStream);
        // exhaustive
        globals.outputEvidence_fileStream << "K,logEvidence_exhaustive";
        // harmonic
        for (int mainRep=0; mainRep<globals.mainRepeats; mainRep++) {
            globals.outputEvidence_fileStream << ",logEvidence_harmonic_rep" << mainRep+1;
        }
        globals.outputEvidence_fileStream << ",logEvidence_harmonic_grandMean,logEvidence_harmonic_grandSE";
        // structure
        for (int mainRep=0; mainRep<globals.mainRepeats; mainRep++) {
            globals.outputEvidence_fileStream << ",logEvidence_structure_rep" << mainRep+1;
        }
        globals.outputEvidence_fileStream << ",logEvidence_structure_grandMean,logEvidence_structure_grandSE";
        // TI
        globals.outputEvidence_fileStream << ",logEvidence_TI,logEvidence_TI_SE";
        
        globals.outputEvidence_fileStream << "\n";
        globals.outputEvidence_fileStream.flush();
        
    }
    
    // ComparisonStatistics
    if (globals.outputComparisonStatistics_on) {
        // open file stream
        globals.outputComparisonStatistics_fileStream = safe_ofstream(globals.outputComparisonStatistics_filePath, globals.outputLog_on, globals.outputLog_fileStream);
        // fixed headers
        globals.outputComparisonStatistics_fileStream << "K,AIC,BIC,DIC_S,DIC_G";
        
        globals.outputComparisonStatistics_fileStream << "\n";
        globals.outputComparisonStatistics_fileStream.flush();
        
    }
    
    // Evanno
    if (globals.outputEvanno_on) {
        // open file stream
        globals.outputEvanno_fileStream = safe_ofstream(globals.outputEvanno_filePath, globals.outputLog_on, globals.outputLog_fileStream);
        // fixed headers
        globals.outputEvanno_fileStream << "K,delta_K";
        
        globals.outputEvanno_fileStream << "\n";
        globals.outputEvanno_fileStream.flush();
        
    }
    
    // EvidenceDetails
    if (globals.outputEvidenceDetails_on) {
        // open file stream
        globals.outputEvidenceDetails_fileStream = safe_ofstream(globals.outputEvidenceDetails_filePath, globals.outputLog_on, globals.outputLog_fileStream);
        // fixed headers
        globals.outputEvidenceDetails_fileStream << "K";
        // structure mean and variance
        for (int mainRep=0; mainRep<globals.mainRepeats; mainRep++) {
            globals.outputEvidenceDetails_fileStream << ",structure_loglike_mean_rep" << mainRep+1;
        }
        for (int mainRep=0; mainRep<globals.mainRepeats; mainRep++) {
            globals.outputEvidenceDetails_fileStream << ",structure_loglike_var_rep" << mainRep+1;
        }
        // TI mean and standard error
        if (globals.thermodynamic_on) {
            for (int TIrep=0; TIrep<globals.thermodynamicRungs; TIrep++) {
                globals.outputEvidenceDetails_fileStream << ",TIpoint_mean_rung" << TIrep+1;
            }
            for (int TIrep=0; TIrep<globals.thermodynamicRungs; TIrep++) {
                globals.outputEvidenceDetails_fileStream << ",TIpoint_SE_rung" << TIrep+1;
            }
        }
        
        globals.outputEvidenceDetails_fileStream << "\n";
        globals.outputEvidenceDetails_fileStream.flush();
        
    }
    
    // PosteriorGrouping
    if (globals.outputPosteriorGrouping_on) {
        // open file stream
        globals.outputPosteriorGrouping_fileStream = safe_ofstream(globals.outputPosteriorGrouping_filePath, globals.outputLog_on, globals.outputLog_fileStream);
        // fixed headers
        globals.outputPosteriorGrouping_fileStream << "K,mainRep,MCMCsample";
        // non-admixture version
        if (!globals.admix_on) {
            for (int i=0; i<globals.n; i++) {
                globals.outputPosteriorGrouping_fileStream << ",ind" << i+1;
            }
        }
        // admixture version
        if (globals.admix_on) {
            for (int i=0; i<globals.n; i++) {
                for (int l=0; l<globals.loci; l++) {
                    for (int p=0; p<globals.ploidy_vec[i]; p++) {
                        globals.outputPosteriorGrouping_fileStream << ",ind" << i+1 << "_loc" << l+1 << "_copy" << p+1;
                    }
                }
            }
        }
        globals.outputPosteriorGrouping_fileStream << "\n";
        globals.outputPosteriorGrouping_fileStream.flush();
    }
    
    // MaxLike_alleleFreqs
    if (globals.outputMaxLike_alleleFreqs_on) {
        // open file stream
        globals.outputMaxLike_alleleFreqs_fileStream = safe_ofstream(globals.outputMaxLike_alleleFreqs_filePath, globals.outputLog_on, globals.outputLog_fileStream);
        // fixed headers
        globals.outputMaxLike_alleleFreqs_fileStream << "K,deme,locus,allele,frequency";
        
        globals.outputMaxLike_alleleFreqs_fileStream << "\n";
        globals.outputMaxLike_alleleFreqs_fileStream.flush();
    }
    
    // MaxLike_admixFreqs
    if (globals.outputMaxLike_admixFreqs_on) {
        // open file stream
        globals.outputMaxLike_admixFreqs_fileStream = safe_ofstream(globals.outputMaxLike_admixFreqs_filePath, globals.outputLog_on, globals.outputLog_fileStream);
        // fixed headers
        globals.outputMaxLike_admixFreqs_fileStream << "K,ind,deme,admixture_proportion";
        
        globals.outputMaxLike_admixFreqs_fileStream << "\n";
        globals.outputMaxLike_admixFreqs_fileStream.flush();
    }
    
}

//------------------------------------------------
// write results to Evidence file
void printEvidence(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    globals.outputEvidence_fileStream << K;
    
    // exhaustive
    globals.outputEvidence_fileStream << "," << process_nan(globals.logEvidence_exhaustive[Kindex]);
    // harmonic
    for (int mainRep=0; mainRep<globals.mainRepeats; mainRep++) {
        globals.outputEvidence_fileStream << "," << process_nan(globals.logEvidence_harmonic[Kindex][mainRep]);
    }
    globals.outputEvidence_fileStream << "," << process_nan(globals.logEvidence_harmonic_grandMean[Kindex]);
    globals.outputEvidence_fileStream << "," << process_nan(globals.logEvidence_harmonic_grandSE[Kindex]);
    // structure
    for (int mainRep=0; mainRep<globals.mainRepeats; mainRep++) {
        globals.outputEvidence_fileStream << "," << process_nan(globals.logEvidence_structure[Kindex][mainRep]);
    }
    globals.outputEvidence_fileStream << "," << process_nan(globals.logEvidence_structure_grandMean[Kindex]);
    globals.outputEvidence_fileStream << "," << process_nan(globals.logEvidence_structure_grandSE[Kindex]);
    // TI
    globals.outputEvidence_fileStream << "," << process_nan(globals.logEvidence_TI[Kindex]);
    globals.outputEvidence_fileStream << "," << process_nan(globals.logEvidence_TI_SE[Kindex]);
    
    globals.outputEvidence_fileStream << "\n";
    globals.outputEvidence_fileStream.flush();
    
}

//------------------------------------------------
// write results to EvidenceDetails file
void printEvidenceDetails(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    globals.outputEvidenceDetails_fileStream << K;
    
    // Structure estimator details
    for (int mainRep=0; mainRep<globals.mainRepeats; mainRep++) {
        globals.outputEvidenceDetails_fileStream << "," << process_nan(globals.structure_loglike_mean[Kindex][mainRep]);
    }
    for (int mainRep=0; mainRep<globals.mainRepeats; mainRep++) {
        globals.outputEvidenceDetails_fileStream << "," << process_nan(globals.structure_loglike_var[Kindex][mainRep]);
    }
    
    // thermodynamic integral estimator details
    if (globals.thermodynamic_on) {
        for (int TIrep=0; TIrep<globals.thermodynamicRungs; TIrep++) {
            globals.outputEvidenceDetails_fileStream << "," << process_nan(globals.TIpoint_mean[Kindex][TIrep]);
        }
        for (int TIrep=0; TIrep<globals.thermodynamicRungs; TIrep++) {
            globals.outputEvidenceDetails_fileStream << "," << process_nan(globals.TIpoint_SE[Kindex][TIrep]);
        }
    }
    
    globals.outputEvidenceDetails_fileStream << "\n";
    globals.outputEvidenceDetails_fileStream.flush();
}

//------------------------------------------------
// write results to EvidenceNormalised file
void printEvidenceNormalised(globals &globals) {
    
    cout << "Calculating normalised results...\n\n";
    
    // normalise exhaustive results
    if (globals.exhaustive_on)
        globals.posterior_exhaustive = normalise_log(globals.logEvidence_exhaustive);
    
    // normalise harmonic mean results
    if (globals.mainRepeats==1) {
        globals.posterior_harmonic_mean = normalise_log(globals.logEvidence_harmonic_grandMean);
    } else {
        normalise_log_sim(globals.posterior_harmonic_mean, globals.posterior_harmonic_LL, globals.posterior_harmonic_UL, globals.logEvidence_harmonic_grandMean, globals.logEvidence_harmonic_grandSE, int(1e6));
    }
    
    // normalise structure estimator results
    if (globals.mainRepeats==1) {
        globals.posterior_structure_mean = normalise_log(globals.logEvidence_structure_grandMean);
    } else {
        normalise_log_sim(globals.posterior_structure_mean, globals.posterior_structure_LL, globals.posterior_structure_UL, globals.logEvidence_structure_grandMean, globals.logEvidence_structure_grandSE, int(1e6));
    }
    
    // normalise thermodynamic integral estimator results
    if (globals.thermodynamic_on)
        normalise_log_sim(globals.posterior_TI_mean, globals.posterior_TI_LL, globals.posterior_TI_UL, globals.logEvidence_TI, globals.logEvidence_TI_SE, int(1e6));
    
    // open file stream
    globals.outputEvidenceNormalised_fileStream = safe_ofstream(globals.outputEvidenceNormalised_filePath, globals.outputLog_on, globals.outputLog_fileStream);
    
    // print headers
    globals.outputEvidenceNormalised_fileStream << "K,posterior_exhaustive,posterior_harmonic_mean,posterior_harmonic_LL,posterior_harmonic_UL,posterior_structure_mean,posterior_structure_LL,posterior_structure_UL,posterior_TI_mean,posterior_TI_LL,posterior_TI_UL\n";
    
    // loop over K
    for (int Kindex=0; Kindex<=(globals.Kmax-globals.Kmin); Kindex++) {
        int K = globals.Kmin+Kindex;
        
        // K
        globals.outputEvidenceNormalised_fileStream << K;
        
        // exhaustive results normalised
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_exhaustive[Kindex]);
        
        // harmonic mean results normalised
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_harmonic_mean[Kindex]);
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_harmonic_LL[Kindex]);
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_harmonic_UL[Kindex]);
        
        // Structure estimator results normalised
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_structure_mean[Kindex]);
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_structure_LL[Kindex]);
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_structure_UL[Kindex]);
        
        // thermodynamic integral estimator results normalised
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_TI_mean[Kindex]);
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_TI_LL[Kindex]);
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_TI_UL[Kindex]);
        
        globals.outputEvidenceNormalised_fileStream << "\n";
        
    }
    globals.outputEvidenceNormalised_fileStream.flush();
}

//------------------------------------------------
// write max-like allele frequencies to file
void printMaxLike_alleleFreqs(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // print values
    for (int k=0; k<K; k++) {
        for (int l=0; l<globals.loci; l++) {
            for (int j=0; j<globals.J[l]; j++) {
                globals.outputMaxLike_alleleFreqs_fileStream << K << "," << k+1 << "," << l+1 << "," << globals.uniqueAlleles[l][j] << "," << process_nan(globals.max_alleleFreqs[k][l][j]) << "\n";
            }
        }
    }
    globals.outputMaxLike_alleleFreqs_fileStream.flush();
    
}

//------------------------------------------------
// write max-like admixture frequencies to file
void printMaxLike_admixFreqs(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // print values
    for (int i=0; i<globals.n; i++) {
        for (int k=0; k<K; k++) {
            globals.outputMaxLike_admixFreqs_fileStream << K << "," << i+1 << "," << k+1 << "," << process_nan(globals.max_admixFreqs[i][k]) << "\n";
        }
    }
    globals.outputMaxLike_alleleFreqs_fileStream.flush();
    
}

//------------------------------------------------
// write comparison statistics to file
void printComparisonStatistics(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    globals.outputComparisonStatistics_fileStream << K << "," << process_nan(globals.AIC[Kindex]) << "," << process_nan(globals.BIC[Kindex]) << "," << process_nan(globals.DIC_Spiegelhalter[Kindex]) << "," << process_nan(globals.DIC_Gelman[Kindex]);

    globals.outputComparisonStatistics_fileStream << "\n";
    globals.outputComparisonStatistics_fileStream.flush();
    
}

//------------------------------------------------
// write Evanno's delta K to file
void printEvanno(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    if (K==1 || K==globals.Kmax) {
        globals.outputEvanno_fileStream << K << ",NA";
    } else {
        globals.outputEvanno_fileStream << K << "," << globals.delta_K[Kindex];
    }
    
    globals.outputEvanno_fileStream << "\n";
    globals.outputEvanno_fileStream.flush();
}

//------------------------------------------------
// write gene-level Qmatrix file (MCMCobject_admixture only)
void printQmatrix_gene(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // open outputQmatrix_gene filestream for this K
    string s = globals.outputQmatrix_gene_filePath;
    s.erase(end(s)-4,end(s));
    s = s + "_K" + to_string((long long)K) + ".csv";
    globals.outputQmatrix_gene_fileStream = safe_ofstream(s,globals.outputLog_on, globals.outputLog_fileStream);
    
    // write headers
    globals.outputQmatrix_gene_fileStream << "index,label";
    if (globals.popCol_on)
        globals.outputQmatrix_gene_fileStream << ",given_population";
    globals.outputQmatrix_gene_fileStream << ",locus,gene_copy";
    for (int k=0; k<K; k++)
        globals.outputQmatrix_gene_fileStream << ",deme" << k+1;
    globals.outputQmatrix_gene_fileStream << "\n";
    
    // write results
    int groupIndex=-1;
    for (int ind=0; ind<globals.n; ind++) {
        for (int l=0; l<globals.loci; l++) {
            for (int p=0; p<globals.ploidy_vec[ind]; p++) {
                groupIndex++;
                
                // individual index
                globals.outputQmatrix_gene_fileStream << ind+1 << ",";
                
                // individual label
                globals.outputQmatrix_gene_fileStream << globals.indLabels_vec[ind] << ",";
                
                // population data
                if (globals.popCol_on)
                    globals.outputQmatrix_gene_fileStream << globals.pop_vec[ind] << ",";
                
                // locus
                globals.outputQmatrix_gene_fileStream << l+1 << ",";
                
                // gene copy
                globals.outputQmatrix_gene_fileStream << p+1 << ",";
                
                // probabilities
                for (int k=0; k<K; k++) {
                    char * buffer = new char[256];
                    sprintf(buffer, "%.3f", round(globals.Qmatrix_gene[Kindex][groupIndex][k]*1000)/1000.0);
                    string s = buffer;
                    delete [] buffer;
                    globals.outputQmatrix_gene_fileStream << s;
                    if (k<(K-1))
                        globals.outputQmatrix_gene_fileStream << ",";
                }
                
                globals.outputQmatrix_gene_fileStream << "\n";
                
            }
        }
    }
    globals.outputQmatrix_gene_fileStream.flush();
}

//------------------------------------------------
// write standard error of gene-level Qmatrix (MCMCobject_admixture only)
void printQmatrixError_gene(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // open outputQmatrixError_gene filestream for this K
    string s = globals.outputQmatrixError_gene_filePath;
    s.erase(end(s)-4,end(s));
    s = s + "_K" + to_string((long long)K) + ".csv";
    globals.outputQmatrixError_gene_fileStream = safe_ofstream(s,globals.outputLog_on, globals.outputLog_fileStream);
    
    // write headers
    globals.outputQmatrixError_gene_fileStream << "index,label";
    if (globals.popCol_on)
        globals.outputQmatrixError_gene_fileStream << ",given_population";
    globals.outputQmatrixError_gene_fileStream << ",locus,gene_copy";
    for (int k=0; k<K; k++)
        globals.outputQmatrixError_gene_fileStream << ",deme" << k+1;
    globals.outputQmatrixError_gene_fileStream << "\n";
    
    // write results
    int groupIndex=-1;
    for (int ind=0; ind<globals.n; ind++) {
        for (int l=0; l<globals.loci; l++) {
            for (int p=0; p<globals.ploidy_vec[ind]; p++) {
                groupIndex++;
                
                // individual index
                globals.outputQmatrixError_gene_fileStream << ind+1 << ",";
                
                // individual label
                globals.outputQmatrixError_gene_fileStream << globals.indLabels_vec[ind] << ",";
                
                // population data
                if (globals.popCol_on)
                    globals.outputQmatrixError_gene_fileStream << globals.pop_vec[ind] << ",";
                
                // locus
                globals.outputQmatrixError_gene_fileStream << l+1 << ",";
                
                // gene copy
                globals.outputQmatrixError_gene_fileStream << p+1 << ",";
                
                // probabilities
                for (int k=0; k<K; k++) {
                    char * buffer = new char[256];
                    sprintf(buffer, "%.3f", round(globals.QmatrixError_gene[Kindex][groupIndex][k]*1000)/1000.0);
                    string s = buffer;
                    delete [] buffer;
                    globals.outputQmatrixError_gene_fileStream << s;
                    if (k<(K-1))
                        globals.outputQmatrixError_gene_fileStream << ",";
                }
                
                globals.outputQmatrixError_gene_fileStream << "\n";
                
            }
        }
    }
    globals.outputQmatrixError_gene_fileStream.flush();
}

//------------------------------------------------
// write individual-level Qmatrix file (overloaded for MCMCobject_noAdmixture and MCMCobject_admixture)
void printQmatrix_ind(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    if (globals.outputQmatrix_structureFormat_on) {
        
        // open outputQmatrix_ind filestream for this K
        string s = globals.outputQmatrix_ind_filePath;
        s.erase(end(s)-4,end(s));
        s = s + "_K" + to_string((long long)K) + ".txt";
        globals.outputQmatrix_ind_fileStream = safe_ofstream(s,globals.outputLog_on, globals.outputLog_fileStream);
        
        // write headers
        /*
        globals.outputQmatrix_ind_fileStream << "        Label (%Miss) ";
        if (globals.popCol_on)
            globals.outputQmatrix_ind_fileStream << "Pop";
        globals.outputQmatrix_ind_fileStream << ":  Inferred clusters\n";
        */
        
        // write results
        for (int i=0; i<globals.n; i++) {
            
            // individual index
            if ((i+1)<100)
                globals.outputQmatrix_ind_fileStream << " ";
            if ((i+1)<10)
                globals.outputQmatrix_ind_fileStream << " ";
            globals.outputQmatrix_ind_fileStream << i+1;
            
            // individual label
            int stringLength = int(globals.indLabels_vec[i].length());
            if (stringLength<8) {
                for (int j=0; j<(8-stringLength); j++) {
                    globals.outputQmatrix_ind_fileStream << " ";
                }
            }
            globals.outputQmatrix_ind_fileStream << " ";
            string s = globals.indLabels_vec[i];
            if (stringLength>11)
                s.erase(begin(s)+11,end(s));
            globals.outputQmatrix_ind_fileStream << s;
            
            // missing data
            double missPercent = round(globals.missing_vec[i]/double(globals.loci)/double(globals.ploidy_vec[i])*100);
            if (missPercent<10)
                globals.outputQmatrix_ind_fileStream << " ";
            globals.outputQmatrix_ind_fileStream << "  (" << missPercent << ")   ";
            
            // population data
            if (globals.popCol_on) {
                if (globals.pop_vec[i].length()==1)
                    globals.outputQmatrix_ind_fileStream << " ";
                globals.outputQmatrix_ind_fileStream << " " << globals.pop_vec[i] << " ";
            }
            
            globals.outputQmatrix_ind_fileStream << ": ";
            
            // probabilities
            for (int k=0; k<K; k++) {
                char * buffer = new char[256];
                sprintf(buffer, "%.3f", round(globals.Qmatrix_ind[Kindex][i][k]*1000)/1000.0);
                string s = buffer;
                delete [] buffer;
                globals.outputQmatrix_ind_fileStream << " " << s;
            }
            globals.outputQmatrix_ind_fileStream << "\n";
            
        }
        
    } // end of if outputQmatrix_structureFormat_on=true
    
    // open outputQmatrix_ind filestream for this K
    string s = globals.outputQmatrix_ind_filePath;
    s.erase(end(s)-4,end(s));
    s = s + "_K" + to_string((long long)K) + ".csv";
    globals.outputQmatrix_ind_fileStream = safe_ofstream(s,globals.outputLog_on, globals.outputLog_fileStream);
    
    // write headers
    globals.outputQmatrix_ind_fileStream << "index,label";
    if (globals.popCol_on)
        globals.outputQmatrix_ind_fileStream << ",given_population";
    for (int k=0; k<K; k++)
        globals.outputQmatrix_ind_fileStream << ",deme" << k+1;
    globals.outputQmatrix_ind_fileStream << "\n";
    
    // write results
    for (int i=0; i<globals.n; i++) {
        
        // individual index
        globals.outputQmatrix_ind_fileStream << i+1;
        globals.outputQmatrix_ind_fileStream << ",";
        
        // individual label
        globals.outputQmatrix_ind_fileStream << globals.indLabels_vec[i];
        globals.outputQmatrix_ind_fileStream << ",";
        
        // population data
        if (globals.popCol_on) {
            globals.outputQmatrix_ind_fileStream << globals.pop_vec[i];
            globals.outputQmatrix_ind_fileStream << ",";
        }
        
        // probabilities
        for (int k=0; k<K; k++) {
            char * buffer = new char[256];
            sprintf(buffer, "%.3f", round(globals.Qmatrix_ind[Kindex][i][k]*1000)/1000.0);
            string s = buffer;
            delete [] buffer;
            globals.outputQmatrix_ind_fileStream << s;
            if (k<(K-1))
                globals.outputQmatrix_ind_fileStream << ",";
        }
        
        globals.outputQmatrix_ind_fileStream << "\n";
        
    }
    
    globals.outputQmatrix_ind_fileStream.flush();
}

//------------------------------------------------
// write standard error of individual-level Qmatrix (overloaded for MCMCobject_noAdmixture and MCMCobject_admixture)
void printQmatrixError_ind(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // open outputQmatrixError_ind filestream for this K
    string s = globals.outputQmatrixError_ind_filePath;
    s.erase(end(s)-4,end(s));
    s = s + "_K" + to_string((long long)K) + ".csv";
    globals.outputQmatrixError_ind_fileStream = safe_ofstream(s,globals.outputLog_on, globals.outputLog_fileStream);
    
    // write headers
    globals.outputQmatrixError_ind_fileStream << "index,label";
    if (globals.popCol_on)
        globals.outputQmatrixError_ind_fileStream << ",given_population";
    for (int k=0; k<K; k++)
        globals.outputQmatrixError_ind_fileStream << ",deme" << k+1;
    globals.outputQmatrixError_ind_fileStream << "\n";
    
    // write results
    for (int i=0; i<globals.n; i++) {
        
        // individual index
        globals.outputQmatrixError_ind_fileStream << i+1;
        globals.outputQmatrixError_ind_fileStream << ",";
        
        // individual label
        globals.outputQmatrixError_ind_fileStream << globals.indLabels_vec[i];
        globals.outputQmatrixError_ind_fileStream << ",";
        
        // population data
        if (globals.popCol_on) {
            globals.outputQmatrixError_ind_fileStream << globals.pop_vec[i];
            globals.outputQmatrixError_ind_fileStream << ",";
        }
        
        // standard errors
        for (int k=0; k<K; k++) {
            char * buffer = new char[256];
            sprintf(buffer, "%.6f", round(globals.QmatrixError_ind[Kindex][i][k]*1e6)/1e6);
            string s = buffer;
            delete [] buffer;
            globals.outputQmatrixError_ind_fileStream << s;
            if (k<(K-1))
                globals.outputQmatrixError_ind_fileStream << ",";
        }
        
        globals.outputQmatrixError_ind_fileStream << "\n";      
    }
    
    globals.outputQmatrixError_ind_fileStream.flush();
}

//------------------------------------------------
// write population-level Qmatrix file (overloaded for MCMCobject_noAdmixture and MCMCobject_admixture)
void printQmatrix_pop(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    if (globals.outputQmatrix_structureFormat_on) {
        
        // open outputQmatrix_pop filestream for this K
        string s = globals.outputQmatrix_pop_filePath;
        s.erase(end(s)-4,end(s));
        s = s + "_K" + to_string((long long)K) + ".txt";
        globals.outputQmatrix_pop_fileStream = safe_ofstream(s,globals.outputLog_on, globals.outputLog_fileStream);
        
        // write headers
        /*
        globals.outputQmatrix_pop_fileStream << "Given    Inferred Clusters       Number of\n";
        globals.outputQmatrix_pop_fileStream << " Pop       1";
        for (int k=1; k<K; k++) {
            if ((k+1)<10)
                globals.outputQmatrix_pop_fileStream << " ";
            globals.outputQmatrix_pop_fileStream << "     " << k+1;
        }
        globals.outputQmatrix_pop_fileStream << "      Individuals\n\n";
        */
        
        // write results
        for (int i=0; i<int(globals.uniquePops.size()); i++) {
            
            // given pop
            if (globals.uniquePops[i].length()<3)
                globals.outputQmatrix_pop_fileStream << " ";
            if (globals.uniquePops[i].length()<2)
                globals.outputQmatrix_pop_fileStream << " ";
            globals.outputQmatrix_pop_fileStream << globals.uniquePops[i] << ":     ";
            
            // probabilities
            for (int k=0; k<K; k++) {
                char * buffer = new char[256];
                sprintf(buffer, "%.3f", round(globals.Qmatrix_pop[Kindex][i][k]*1000)/1000.0);
                string s = buffer;
                delete [] buffer;
                globals.outputQmatrix_pop_fileStream << s << " ";
            }
            
            // number of individuals in each population
            globals.outputQmatrix_pop_fileStream << "     ";
            if (globals.uniquePop_counts[i]<100)
                globals.outputQmatrix_pop_fileStream << " ";
            if (globals.uniquePop_counts[i]<10)
                globals.outputQmatrix_pop_fileStream << " ";
            globals.outputQmatrix_pop_fileStream << globals.uniquePop_counts[i];
            
            globals.outputQmatrix_pop_fileStream << "\n";
        }
        
    } // end of if outputQmatrix_structureFormat_on=true
        
    // open outputQmatrix_pop filestream for this K
    string s = globals.outputQmatrix_pop_filePath;
    s.erase(end(s)-4,end(s));
    s = s + "_K" + to_string((long long)K) + ".csv";
    globals.outputQmatrix_pop_fileStream = safe_ofstream(s,globals.outputLog_on, globals.outputLog_fileStream);
    
    // write headers
    globals.outputQmatrix_pop_fileStream << "given_population,members";
    for (int k=0; k<K; k++)
        globals.outputQmatrix_pop_fileStream << ",deme" << k+1;
    globals.outputQmatrix_pop_fileStream << "\n";
    
    // write results
    for (int i=0; i<int(globals.uniquePops.size()); i++) {
        
        // given population and number of individuals
        globals.outputQmatrix_pop_fileStream << globals.uniquePops[i] << "," << globals.uniquePop_counts[i];
        
        // probabilities
        for (int k=0; k<K; k++) {
            char * buffer = new char[256];
            sprintf(buffer, "%.3f", round(globals.Qmatrix_pop[Kindex][i][k]*1000)/1000.0);
            string s = buffer;
            delete [] buffer;
            globals.outputQmatrix_pop_fileStream << "," << s;
        }
        
        globals.outputQmatrix_pop_fileStream << "\n";
        
    }
    
    globals.outputQmatrix_pop_fileStream.flush();
}

//------------------------------------------------
// write standard error of population-level Qmatrix (overloaded for MCMCobject_noAdmixture and MCMCobject_admixture)
void printQmatrixError_pop(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    // open outputQmatrixError_pop filestream for this K
    string s = globals.outputQmatrixError_pop_filePath;
    s.erase(end(s)-4,end(s));
    s = s + "_K" + to_string((long long)K) + ".csv";
    globals.outputQmatrixError_pop_fileStream = safe_ofstream(s,globals.outputLog_on, globals.outputLog_fileStream);
    
    // write headers
    globals.outputQmatrixError_pop_fileStream << "given_population,individuals";
    for (int k=0; k<K; k++)
        globals.outputQmatrixError_pop_fileStream << ",deme" << k+1;
    globals.outputQmatrixError_pop_fileStream << "\n";
    
    // write results
    for (int i=0; i<int(globals.uniquePops.size()); i++) {
        
        // given population and number of individuals
        globals.outputQmatrixError_pop_fileStream << globals.uniquePops[i] << "," << globals.uniquePop_counts[i];
        
        // standard errors
        for (int k=0; k<K; k++) {
            char * buffer = new char[256];
            sprintf(buffer, "%.6f", round(globals.QmatrixError_pop[Kindex][i][k]*1e6)/1e6);
            string s = buffer;
            delete [] buffer;
            globals.outputQmatrixError_pop_fileStream << "," << s;
        }
        
        globals.outputQmatrixError_pop_fileStream << "\n";
        
    }
    
    globals.outputQmatrixError_pop_fileStream.flush();
    
}
