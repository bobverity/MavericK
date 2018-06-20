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
        // TI
        globals.outputEvidence_fileStream << ",logEvidence_TI,logEvidence_TI_SE";
        
        globals.outputEvidence_fileStream << "\n";
        globals.outputEvidence_fileStream.flush();
        
    }
    
    // EvidenceDetails
    if (globals.outputEvidenceDetails_on) {
        // open file stream
        globals.outputEvidenceDetails_fileStream = safe_ofstream(globals.outputEvidenceDetails_filePath, globals.outputLog_on, globals.outputLog_fileStream);
        // fixed headers
        globals.outputEvidenceDetails_fileStream << "K";
        // TI mean and standard error
        for (int TIrep=0; TIrep<globals.rungs; TIrep++) {
            globals.outputEvidenceDetails_fileStream << ",TIpoint_mean_rung" << TIrep+1;
        }
        for (int TIrep=0; TIrep<globals.rungs; TIrep++) {
            globals.outputEvidenceDetails_fileStream << ",TIpoint_SE_rung" << TIrep+1;
        }
        globals.outputEvidenceDetails_fileStream << "\n";
        globals.outputEvidenceDetails_fileStream.flush();
        
    }
    
}

//------------------------------------------------
// write results to Evidence file
void printEvidence(globals &globals, int Kindex) {
    int K = globals.Kmin+Kindex;
    
    globals.outputEvidence_fileStream << K;
    
    // exhaustive
    globals.outputEvidence_fileStream << "," << process_nan(globals.logEvidence_exhaustive[Kindex]);
    
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
    
    // thermodynamic integral estimator details
    for (int TIrep=0; TIrep<globals.rungs; TIrep++) {
        globals.outputEvidenceDetails_fileStream << "," << process_nan(globals.TIpoint_mean[Kindex][TIrep]);
    }
    for (int TIrep=0; TIrep<globals.rungs; TIrep++) {
        globals.outputEvidenceDetails_fileStream << "," << process_nan(globals.TIpoint_SE[Kindex][TIrep]);
    }
    
    globals.outputEvidenceDetails_fileStream << "\n";
    globals.outputEvidenceDetails_fileStream.flush();
}

//------------------------------------------------
// write results to EvidenceNormalised file
void printEvidenceNormalised(globals &globals) {
    
    cout << "Calculating normalised results...\n\n";
    
    // normalise exhaustive results
    if (globals.exhaustive_on) {
        globals.posterior_exhaustive = normalise_log(globals.logEvidence_exhaustive);
    }
    
    // normalise thermodynamic integral estimator results
    normalise_log_sim(globals.posterior_TI_mean, globals.posterior_TI_LL, globals.posterior_TI_UL, globals.logEvidence_TI, globals.logEvidence_TI_SE, int(1e6));
    
    // open file stream
    globals.outputEvidenceNormalised_fileStream = safe_ofstream(globals.outputEvidenceNormalised_filePath, globals.outputLog_on, globals.outputLog_fileStream);
    
    // print headers
    globals.outputEvidenceNormalised_fileStream << "K,posterior_exhaustive,posterior_TI_mean,posterior_TI_LL,posterior_TI_UL\n";
    
    // loop over K
    for (int Kindex=0; Kindex<=(globals.Kmax-globals.Kmin); Kindex++) {
        int K = globals.Kmin+Kindex;
        
        // K
        globals.outputEvidenceNormalised_fileStream << K;
        
        // exhaustive results normalised
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_exhaustive[Kindex]);
        
        // thermodynamic integral estimator results normalised
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_TI_mean[Kindex]);
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_TI_LL[Kindex]);
        globals.outputEvidenceNormalised_fileStream << "," << process_nan(globals.posterior_TI_UL[Kindex]);
        
        globals.outputEvidenceNormalised_fileStream << "\n";
        
    }
    globals.outputEvidenceNormalised_fileStream.flush();
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
