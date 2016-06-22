//
//  MavericK
//  readIn.cpp
//
//  Created: Bob on 22/09/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "readIn.h"

using namespace std;

//------------------------------------------------
// write to file if writeToFile is true
void writeToFile(string message, bool writeToFile, ofstream &fileStream) {
    if (writeToFile) {
        fileStream << message;
        fileStream.flush();
    }
}

//------------------------------------------------
// write message to screen, and to file if writeToFile is true
void coutAndLog(string message, bool writeToFile, ofstream &fileStream) {
    cout << message;
    if (writeToFile) {
        fileStream << message;
        fileStream.flush();
    }
}

//------------------------------------------------
// write error to screen, and to file if writeToFile is true
void cerrAndLog(string message, bool writeToFile, ofstream &logFileStream) {
    cerr << message;
    if (writeToFile) {
        logFileStream << message;
        logFileStream.flush();
    }
}

//------------------------------------------------
// safely open input file stream, otherwise error. Option to write errors to file if writeErrorToFile is true.
ifstream safe_ifstream(string fileName, bool writeToFile, ofstream &logFileStream) {
    ifstream fileStream(fileName);
    if (!fileStream.is_open()) {
        cerrAndLog("\nError: failed to read from file: "+fileName+string("\n"), writeToFile, logFileStream);
        exit(1);
    }
    return(fileStream);
}

//------------------------------------------------
// convert string to double if possible, otherwise error
double stringToDouble(string s, bool writeToFile, ofstream &logFileStream) {
    double output;
    istringstream ss(s);
    ss >> output;
    if (!output) {
        cerrAndLog("\nError: cannot convert data element '"+s+string("' to type double\n"), writeToFile, logFileStream);
        exit(1);
    }
    return(output);
}

//------------------------------------------------
// convert double to int (if no remainder), otherwise error
int doubleToInt(double x, bool writeToFile, ofstream &logFileStream) {
    double intpart;
    if (modf(x, &intpart)!=0.0) {
        cerrAndLog("\nError: cannot convert data element "+to_string((long long)x)+string(" to type integer\n"), writeToFile, logFileStream);
        exit(1);
    }
    return(int(intpart));
}

//------------------------------------------------
// equivalent to getline() function, but searches for '\r', '\n' or '\r\n' (therefore Max, PC, Linux compatible)
istream& safe_getline(istream& is, string& t) {
	t.clear();
	istream::sentry se(is, true);
	streambuf* sb = is.rdbuf();
	for (;;) {
		int c = sb->sbumpc();
		switch(c) {
            case '\n':
                return is;
            case '\r':
                if (sb->sgetc()=='\n')
                    sb->sbumpc();
                return is;
            case EOF:
                if (t.empty())
                    is.setstate(ios::eofbit);
                return is;
            default:
                t += (char)c;
		}
	}
}

//------------------------------------------------
// if argument argv[i] is equal to argName then read item [i+1] into savePath
void readPath(string argName, string &savePath, int argc, const char * argv[], int i) {
    if (string(argv[i])==argName) {
        if (i+1<argc) {
            savePath = argv[i+1];
        } else {
            cerr << "Error: there was no argument to the "+argName+string(" option\n");
            exit(1);
        }
    }
}

//------------------------------------------------
// read in parameters as strings
void readParameters(globals &globals) {
    
    cout << "Loading parameters...\n";
    
    // open file stream
    globals.parameters_fileStream = safe_ifstream(globals.parameters_filePath, false, globals.outputLog_fileStream);
    
    // read in parameters file as strings
    vector<string> params;
    string line1;
    while (safe_getline(globals.parameters_fileStream, line1)) {
        istringstream ss(line1);
        string line2;
        while (getline(ss, line2, '\t')) {
            if (line2.size()>0) {
                string line3 = istringstream(line2).str();
                params.push_back(line3);
            }
        }
    }
    
    // close file stream
    globals.parameters_fileStream.close();
    
    // extract paths from file
    for (int i=1; i<int(params.size()); i++) {
        
        if (params[i]=="outputRoot" && i+1<int(params.size()))
            globals.outputRoot_fileName = params[i+1];
        
        if (params[i]=="data" && i+1<int(params.size()))
            globals.data_fileName = params[i+1];
        
        if (params[i]=="outputLog" && i+1<int(params.size()))
            globals.outputLog_fileName = params[i+1];
        
        if (params[i]=="outputLikelihood" && i+1<int(params.size()))
            globals.outputLikelihood_fileName = params[i+1];
        
        if (params[i]=="outputQmatrix_ind" && i+1<int(params.size()))
            globals.outputQmatrix_ind_fileName = params[i+1];
        
        if (params[i]=="outputQmatrix_pop" && i+1<int(params.size()))
            globals.outputQmatrix_pop_fileName = params[i+1];
        
        if (params[i]=="outputQmatrix_gene" && i+1<int(params.size()))
            globals.outputQmatrix_gene_fileName = params[i+1];
        
        if (params[i]=="outputQmatrixError_ind" && i+1<int(params.size()))
            globals.outputQmatrixError_ind_fileName = params[i+1];
        
        if (params[i]=="outputQmatrixError_pop" && i+1<int(params.size()))
            globals.outputQmatrixError_pop_fileName = params[i+1];
        
        if (params[i]=="outputQmatrixError_gene" && i+1<int(params.size()))
            globals.outputQmatrixError_gene_fileName = params[i+1];
        
        if (params[i]=="outputEvidence" && i+1<int(params.size()))
            globals.outputEvidence_fileName = params[i+1];
        
        if (params[i]=="outputEvidenceNormalised" && i+1<int(params.size()))
            globals.outputEvidenceNormalised_fileName = params[i+1];
        
        if (params[i]=="outputEvidenceDetails" && i+1<int(params.size()))
            globals.outputEvidenceDetails_fileName = params[i+1];
        
        if (params[i]=="outputPosteriorGrouping" && i+1<int(params.size()))
            globals.outputPosteriorGrouping_fileName = params[i+1];
        
        if (params[i]=="outputComparisonStatistics" && i+1<int(params.size()))
            globals.outputComparisonStatistics_fileName = params[i+1];
        
        if (params[i]=="outputEvanno" && i+1<int(params.size()))
            globals.outputEvanno_fileName = params[i+1];
        
        if (params[i]=="outputMaxLike_alleleFreqs" && i+1<int(params.size()))
            globals.outputMaxLike_alleleFreqs_fileName = params[i+1];
        
        if (params[i]=="outputMaxLike_admixFreqs" && i+1<int(params.size()))
            globals.outputMaxLike_admixFreqs_fileName = params[i+1];
        
    }
    
    // extract parameter values (as strings) from file
    for (int i=1; i<int(params.size()); i++) {
        
        if (params[i]=="headerRow_on" && i+1<int(params.size()))
            globals.parameterStrings["headerRow_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="popCol_on" && i+1<int(params.size()))
            globals.parameterStrings["popCol_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="ploidyCol_on" && i+1<int(params.size()))
            globals.parameterStrings["ploidyCol_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="ploidy" && i+1<int(params.size()))
            globals.parameterStrings["ploidy"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="missingData" && i+1<int(params.size()))
            globals.parameterStrings["missingData"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="Kmin" && i+1<int(params.size()))
            globals.parameterStrings["Kmin"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="Kmax" && i+1<int(params.size()))
            globals.parameterStrings["Kmax"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="admix_on" && i+1<int(params.size()))
            globals.parameterStrings["admix_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="fixAlpha_on" && i+1<int(params.size()))
            globals.parameterStrings["fixAlpha_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="alpha" && i+1<int(params.size()))
            globals.parameterStrings["alpha"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="alphaPropSD" && i+1<int(params.size()))
            globals.parameterStrings["alphaPropSD"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="exhaustive_on" && i+1<int(params.size()))
            globals.parameterStrings["exhaustive_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="mainRepeats" && i+1<int(params.size()))
            globals.parameterStrings["mainRepeats"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="mainBurnin" && i+1<int(params.size()))
            globals.parameterStrings["mainBurnin"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="mainSamples" && i+1<int(params.size()))
            globals.parameterStrings["mainSamples"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="mainThinning" && i+1<int(params.size()))
            globals.parameterStrings["mainThinning"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="thermodynamic_on" && i+1<int(params.size()))
            globals.parameterStrings["thermodynamic_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="thermodynamicRungs" && i+1<int(params.size()))
            globals.parameterStrings["thermodynamicRungs"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="thermodynamicBurnin" && i+1<int(params.size()))
            globals.parameterStrings["thermodynamicBurnin"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="thermodynamicSamples" && i+1<int(params.size()))
            globals.parameterStrings["thermodynamicSamples"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="thermodynamicThinning" && i+1<int(params.size()))
            globals.parameterStrings["thermodynamicThinning"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="EMalgorithm_on" && i+1<int(params.size()))
            globals.parameterStrings["EMalgorithm_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="EMrepeats" && i+1<int(params.size()))
            globals.parameterStrings["EMrepeats"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="EMiterations" && i+1<int(params.size()))
            globals.parameterStrings["EMiterations"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputLog_on" && i+1<int(params.size()))
            globals.parameterStrings["outputLog_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputLikelihood_on" && i+1<int(params.size()))
            globals.parameterStrings["outputLikelihood_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputQmatrix_ind_on" && i+1<int(params.size()))
            globals.parameterStrings["outputQmatrix_ind_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputQmatrix_pop_on" && i+1<int(params.size()))
            globals.parameterStrings["outputQmatrix_pop_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputQmatrix_gene_on" && i+1<int(params.size()))
            globals.parameterStrings["outputQmatrix_gene_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputQmatrixError_ind_on" && i+1<int(params.size()))
            globals.parameterStrings["outputQmatrixError_ind_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputQmatrixError_pop_on" && i+1<int(params.size()))
            globals.parameterStrings["outputQmatrixError_pop_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputQmatrixError_gene_on" && i+1<int(params.size()))
            globals.parameterStrings["outputQmatrixError_gene_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputEvidence_on" && i+1<int(params.size()))
            globals.parameterStrings["outputEvidence_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputEvidenceNormalised_on" && i+1<int(params.size()))
            globals.parameterStrings["outputEvidenceNormalised_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputEvidenceDetails_on" && i+1<int(params.size()))
            globals.parameterStrings["outputEvidenceDetails_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputPosteriorGrouping_on" && i+1<int(params.size()))
            globals.parameterStrings["outputPosteriorGrouping_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputComparisonStatistics_on" && i+1<int(params.size()))
            globals.parameterStrings["outputComparisonStatistics_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputEvanno_on" && i+1<int(params.size()))
            globals.parameterStrings["outputEvanno_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputMaxLike_alleleFreqs_on" && i+1<int(params.size()))
            globals.parameterStrings["outputMaxLike_alleleFreqs_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputMaxLike_admixFreqs_on" && i+1<int(params.size()))
            globals.parameterStrings["outputMaxLike_admixFreqs_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="outputQmatrix_structureFormat_on" && i+1<int(params.size()))
            globals.parameterStrings["outputQmatrix_structureFormat_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="suppressWarning1_on" && i+1<int(params.size()))
            globals.parameterStrings["suppressWarning1_on"] = pair<string,int>(params[i+1],1);
        
        if (params[i]=="fixLabels_on" && i+1<int(params.size()))
            globals.parameterStrings["fixLabels_on"] = pair<string,int>(params[i+1],1);
        
    }
    
}

//------------------------------------------------
// if argument argv[i] is equal to argName then read item [i+1] into argValue
void readArgument(string argName, globals &globals, int argc, const char * argv[], int i) {
    if (string(argv[i])=="-"+argName) {
        if (i+1<argc) {
            globals.parameterStrings[argName] = pair<string,int>(argv[i+1],2);
        } else {
            cerr << "Error: there was no argument to the "+argName+string(" option\n");
            exit(1);
        }
    }
}

//------------------------------------------------
// read in command line arguments as strings
void readCommandLine(globals &globals, int argc, const char * argv[]) {
    
    // replace file names and paths with user-defined arguments
    for (int i=1; i<argc; i++) {
        readPath("-outputRoot", globals.outputRoot_fileName, argc, argv, i);
        readPath("-data", globals.data_fileName, argc, argv, i);
        readPath("-outputLog", globals.outputLog_fileName, argc, argv, i);
        readPath("-outputLikelihood", globals.outputLikelihood_fileName, argc, argv, i);
        readPath("-outputQmatrix_ind", globals.outputQmatrix_ind_fileName, argc, argv, i);
        readPath("-outputQmatrix_pop", globals.outputQmatrix_pop_fileName, argc, argv, i);
        readPath("-outputQmatrix_gene", globals.outputQmatrix_gene_fileName, argc, argv, i);
        readPath("-outputQmatrixError_ind", globals.outputQmatrixError_ind_fileName, argc, argv, i);
        readPath("-outputQmatrixError_pop", globals.outputQmatrixError_pop_fileName, argc, argv, i);
        readPath("-outputQmatrixError_gene", globals.outputQmatrixError_gene_fileName, argc, argv, i);
        readPath("-outputEvidence", globals.outputEvidence_fileName, argc, argv, i);
        readPath("-outputEvidenceNormalised", globals.outputEvidenceNormalised_fileName, argc, argv, i);
        readPath("-outputEvidenceDetails", globals.outputEvidenceDetails_fileName, argc, argv, i);
        readPath("-outputPosteriorGrouping", globals.outputPosteriorGrouping_fileName, argc, argv, i);
        readPath("-outputComparisonStatistics", globals.outputComparisonStatistics_fileName, argc, argv, i);
        readPath("-outputEvanno", globals.outputEvanno_fileName, argc, argv, i);
        readPath("-outputMaxLike_alleleFreqs", globals.outputMaxLike_alleleFreqs_fileName, argc, argv, i);
        readPath("-outputMaxLike_admixFreqs", globals.outputMaxLike_admixFreqs_fileName, argc, argv, i);
    }
    
    // set file paths
    globals.outputRoot_filePath = globals.masterRoot_filePath + globals.outputRoot_fileName;
    globals.data_filePath = globals.inputRoot_filePath + globals.data_fileName;
    globals.outputLog_filePath = globals.outputRoot_filePath + globals.outputLog_fileName;
    globals.outputLikelihood_filePath = globals.outputRoot_filePath + globals.outputLikelihood_fileName;
    globals.outputQmatrix_ind_filePath = globals.outputRoot_filePath + globals.outputQmatrix_ind_fileName;
    globals.outputQmatrix_pop_filePath = globals.outputRoot_filePath + globals.outputQmatrix_pop_fileName;
    globals.outputQmatrix_gene_filePath = globals.outputRoot_filePath + globals.outputQmatrix_gene_fileName;
    globals.outputQmatrixError_ind_filePath = globals.outputRoot_filePath + globals.outputQmatrixError_ind_fileName;
    globals.outputQmatrixError_pop_filePath = globals.outputRoot_filePath + globals.outputQmatrixError_pop_fileName;
    globals.outputQmatrixError_gene_filePath = globals.outputRoot_filePath + globals.outputQmatrixError_gene_fileName;
    globals.outputEvidence_filePath = globals.outputRoot_filePath + globals.outputEvidence_fileName;
    globals.outputEvidenceNormalised_filePath = globals.outputRoot_filePath + globals.outputEvidenceNormalised_fileName;
    globals.outputEvidenceDetails_filePath = globals.outputRoot_filePath + globals.outputEvidenceDetails_fileName;
    globals.outputPosteriorGrouping_filePath = globals.outputRoot_filePath + globals.outputPosteriorGrouping_fileName;
    globals.outputComparisonStatistics_filePath = globals.outputRoot_filePath + globals.outputComparisonStatistics_fileName;
    globals.outputEvanno_filePath = globals.outputRoot_filePath + globals.outputEvanno_fileName;
    globals.outputMaxLike_alleleFreqs_filePath = globals.outputRoot_filePath + globals.outputMaxLike_alleleFreqs_fileName;
    globals.outputMaxLike_admixFreqs_filePath = globals.outputRoot_filePath + globals.outputMaxLike_admixFreqs_fileName;
    
    // replace parameter values (as strings) with user-defined arguments
    for (int i=1; i<argc; i++) {
        readArgument("headerRow_on", globals, argc, argv, i);
        readArgument("popCol_on", globals, argc, argv, i);
        readArgument("ploidyCol_on", globals, argc, argv, i);
        readArgument("ploidy", globals, argc, argv, i);
        readArgument("missingData", globals, argc, argv, i);
        readArgument("Kmin", globals, argc, argv, i);
        readArgument("Kmax", globals, argc, argv, i);
        readArgument("admix_on", globals, argc, argv, i);
        readArgument("fixAlpha_on", globals, argc, argv, i);
        readArgument("alpha", globals, argc, argv, i);
        readArgument("alphaPropSD", globals, argc, argv, i);
        readArgument("exhaustive_on", globals, argc, argv, i);
        readArgument("mainRepeats", globals, argc, argv, i);
        readArgument("mainBurnin", globals, argc, argv, i);
        readArgument("mainSamples", globals, argc, argv, i);
        readArgument("mainThinning", globals, argc, argv, i);
        readArgument("thermodynamic_on", globals, argc, argv, i);
        readArgument("thermodynamicRungs", globals, argc, argv, i);
        readArgument("thermodynamicBurnin", globals, argc, argv, i);
        readArgument("thermodynamicSamples", globals, argc, argv, i);
        readArgument("thermodynamicThinning", globals, argc, argv, i);
        readArgument("EMalgorithm_on", globals, argc, argv, i);
        readArgument("EMrepeats", globals, argc, argv, i);
        readArgument("EMiterations", globals, argc, argv, i);
        readArgument("outputLog_on", globals, argc, argv, i);
        readArgument("outputLikelihood_on", globals, argc, argv, i);
        readArgument("outputQmatrix_ind_on", globals, argc, argv, i);
        readArgument("outputQmatrix_pop_on", globals, argc, argv, i);
        readArgument("outputQmatrix_gene_on", globals, argc, argv, i);
        readArgument("outputQmatrixError_ind_on", globals, argc, argv, i);
        readArgument("outputQmatrixError_pop_on", globals, argc, argv, i);
        readArgument("outputQmatrixError_gene_on", globals, argc, argv, i);
        readArgument("outputAdmixture_on", globals, argc, argv, i);
        readArgument("outputEvidence_on", globals, argc, argv, i);
        readArgument("outputEvidenceNormalised_on", globals, argc, argv, i);
        readArgument("outputEvidenceDetails_on", globals, argc, argv, i);
        readArgument("outputPosteriorGrouping_on", globals, argc, argv, i);
        readArgument("outputComparisonStatistics_on", globals, argc, argv, i);
        readArgument("outputEvanno_on", globals, argc, argv, i);
        readArgument("outputMaxLike_alleleFreqs_on", globals, argc, argv, i);
        readArgument("outputMaxLike_admixFreqs_on", globals, argc, argv, i);
        readArgument("outputQmatrix_structureFormat_on", globals, argc, argv, i);
        readArgument("suppressWarning1_on", globals, argc, argv, i);
        readArgument("fixLabels_on", globals, argc, argv, i);
    }
    
}

//------------------------------------------------
// check parameters, stored in pair<string,int> format. The input i filters the second element in this pair.
void checkParameters(globals &globals, int i) {
    
    // loop through all elements of paired list
    int defined = 0;
    for (auto it = globals.parameterStrings.begin(); it != globals.parameterStrings.end(); ++it) {
        if (it->second.second==i) {
            defined ++;
            
            if (it->first=="headerRow_on") {
                writeToFile("  headerRow_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.headerRow_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="popCol_on") {
                writeToFile("  popCol_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.popCol_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="ploidyCol_on") {
                writeToFile("  ploidyCol_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.ploidyCol_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="ploidy") {
                writeToFile("  ploidy = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than 0
                checkInteger(it->second.first, globals.ploidy, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                checkGrZero(it->first, globals.ploidy, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="missingData") {
                writeToFile("  missingData = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                globals.missingData = it->second.first;
            }
            if (it->first=="Kmin") {
                writeToFile("  Kmin = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than 0 and less than or equal to Kmax
                checkInteger(it->second.first, globals.Kmin, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                checkGrZero(it->first, globals.Kmin, globals.outputLog_on, globals.outputLog_fileStream);
                if (globals.Kmin>globals.Kmax) {
                    cerrAndLog("\nError: 'Kmin' parameter must be less than or equal to 'Kmax' parameter\n", globals.outputLog_on, globals.outputLog_fileStream);
                    exit(1);
                }
            }
            if (it->first=="Kmax") {
                writeToFile("  Kmax = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than or equal to Kmin
                checkInteger(it->second.first, globals.Kmax, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                if (globals.Kmax<globals.Kmin) {
                    cerrAndLog("\nError: 'Kmax' parameter must be greater than or equal to 'Kmin' parameter\n", globals.outputLog_on, globals.outputLog_fileStream);
                    exit(1);
                }
            }
            if (it->first=="admix_on") {
                writeToFile("  admix_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.admix_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="fixAlpha_on") {
                writeToFile("  fixAlpha_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.fixAlpha_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="alpha") {
                writeToFile("  alpha = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // read in potentially multiple values of alpha
                vector<double> alpha;
                istringstream ss(it->second.first);
                string line1;
                int i=0;
                while (getline(ss, line1, ',')) {
                    alpha.push_back(0);
                    istringstream(line1) >> alpha[i];
                    checkGrZero(it->first, alpha[i], globals.outputLog_on, globals.outputLog_fileStream);
                    i++;
                }
                
                // save vector of values to globals.alpha
                globals.alpha = alpha;
            }
            if (it->first=="alphaPropSD") {
                writeToFile("  alphaPropSD = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // read in potentially multiple values of alphaPropSD
                vector<double> alphaPropSD;
                istringstream ss(it->second.first);
                string line1;
                int i=0;
                while (getline(ss, line1, ',')) {
                    alphaPropSD.push_back(0);
                    istringstream(line1) >> alphaPropSD[i];
                    checkGrZero(it->first, alphaPropSD[i], globals.outputLog_on, globals.outputLog_fileStream);
                    i++;
                }
                
                // save vector of values to globals.alphaPropSD
                globals.alphaPropSD = alphaPropSD;
            }
            if (it->first=="exhaustive_on") {
                writeToFile("  exhaustive_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.exhaustive_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="mainRepeats") {
                writeToFile("  mainRepeats = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than 0
                checkInteger(it->second.first, globals.mainRepeats, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                checkGrZero(it->first, globals.mainRepeats, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="mainBurnin") {
                writeToFile("  mainBurnin = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than or equal to 0
                checkInteger(it->second.first, globals.mainBurnin, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                checkGrEqZero(it->first, globals.mainBurnin, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="mainSamples") {
                writeToFile("  mainSamples = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than 0
                checkInteger(it->second.first, globals.mainSamples, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                checkGrZero(it->first, globals.mainSamples, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="mainThinning") {
                writeToFile("  mainThinning = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than 0
                checkInteger(it->second.first, globals.mainThinning, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                checkGrZero(it->first, globals.mainThinning, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="thermodynamic_on") {
                writeToFile("  thermodynamic_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.thermodynamic_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="thermodynamicRungs") {
                writeToFile("  thermodynamicRungs = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than or equal to 2
                checkInteger(it->second.first, globals.thermodynamicRungs, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                if (globals.thermodynamicRungs<2) {
                    cerrAndLog("\nError: 'thermodynamicRungs' parameter must be greater than or equal to 2\n", globals.outputLog_on, globals.outputLog_fileStream);
                    exit(1);
                }
            }
            if (it->first=="thermodynamicBurnin") {
                writeToFile("  thermodynamicBurnin = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than or equal to 0
                checkInteger(it->second.first, globals.thermodynamicBurnin, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                checkGrEqZero(it->first, globals.thermodynamicBurnin, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="thermodynamicSamples") {
                writeToFile("  thermodynamicSamples = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than 0
                checkInteger(it->second.first, globals.thermodynamicSamples, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                checkGrZero(it->first, globals.thermodynamicSamples, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="thermodynamicThinning") {
                writeToFile("  thermodynamicThinning = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than 0
                checkInteger(it->second.first, globals.thermodynamicThinning, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                checkGrZero(it->first, globals.thermodynamicThinning, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="EMalgorithm_on") {
                writeToFile("  EMalgorithm_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.EMalgorithm_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="EMrepeats") {
                writeToFile("  EMrepeats = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than 0
                checkInteger(it->second.first, globals.EMrepeats, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                checkGrZero(it->first, globals.EMrepeats, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="EMiterations") {
                writeToFile("  EMiterations = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                
                // check that integer greater than 0
                checkInteger(it->second.first, globals.EMiterations, it->first, globals.outputLog_on, globals.outputLog_fileStream);
                checkGrZero(it->first, globals.EMiterations, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputLog_on") {
                writeToFile("  outputLog_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputLog_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputLikelihood_on") {
                writeToFile("  outputLikelihood_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputLikelihood_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputQmatrix_ind_on") {
                writeToFile("  outputQmatrix_ind_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputQmatrix_ind_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputQmatrix_pop_on") {
                writeToFile("  outputQmatrix_pop_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputQmatrix_pop_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputQmatrix_gene_on") {
                writeToFile("  outputQmatrix_gene_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputQmatrix_gene_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputQmatrixError_ind_on") {
                writeToFile("  outputQmatrixError_ind_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputQmatrixError_ind_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputQmatrixError_pop_on") {
                writeToFile("  outputQmatrixError_pop_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputQmatrixError_pop_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputQmatrixError_gene_on") {
                writeToFile("  outputQmatrixError_gene_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputQmatrixError_gene_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputAdmixture_on") {
                writeToFile("  outputAdmixture_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputAdmixture_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputEvidence_on") {
                writeToFile("  outputEvidence_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputEvidence_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputEvidenceNormalised_on") {
                writeToFile("  outputEvidenceNormalised_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputEvidenceNormalised_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputEvidenceDetails_on") {
                writeToFile("  outputEvidenceDetails_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputEvidenceDetails_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputPosteriorGrouping_on") {
                writeToFile("  outputPosteriorGrouping_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputPosteriorGrouping_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputComparisonStatistics_on") {
                writeToFile("  outputComparisonStatistics_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputComparisonStatistics_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputEvanno_on") {
                writeToFile("  outputEvanno_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputEvanno_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputMaxLike_alleleFreqs_on") {
                writeToFile("  outputMaxLike_alleleFreqs_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputMaxLike_alleleFreqs_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputMaxLike_admixFreqs_on") {
                writeToFile("  outputMaxLike_admixFreqs_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputMaxLike_admixFreqs_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="outputQmatrix_structureFormat_on") {
                writeToFile("  outputQmatrix_structureFormat_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.outputQmatrix_structureFormat_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="suppressWarning1_on") {
                writeToFile("  suppressWarning1_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.suppressWarning1_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
            if (it->first=="fixLabels_on") {
                writeToFile("  fixLabels_on = "+it->second.first+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
                checkBoolean(it->second.first, globals.fixLabels_on, it->first, globals.outputLog_on, globals.outputLog_fileStream);
            }
        }
    }
	
    // if no parameters defined by user
    if (defined==0)
        writeToFile("  (none)\n", globals.outputLog_on, globals.outputLog_fileStream);
    
    // print to screen number of parameters of each type
    if (i==0)
        cout << "  " << defined << " parameters set to default values\n";
    else if (i==1)
        cout << "  " << defined << " parameters read in from file\n";
    else if (i==2)
        cout << "  " << defined << " parameters defined on command line\n";
}

//------------------------------------------------
// read in data file
void readData(globals &globals) {
    
    cout << "Loading data file...\n";
    writeToFile("Data properties\n", globals.outputLog_on, globals.outputLog_fileStream);
    
    // open file stream
    globals.data_fileStream = safe_ifstream(globals.data_filePath, globals.outputLog_on, globals.outputLog_fileStream);
    
    // read in raw data file as matrix of strings. Skip header line if necessary
    vector< vector<string> > rawEntries;
    vector<string> rawEntries_line;
    string line1;
    int row=0;
    while (safe_getline(globals.data_fileStream, line1))
    {
        if (line1.size()==0)
            continue;
        row++;
        if (row==1 && globals.headerRow_on)
            continue;
        rawEntries_line.clear();
        istringstream ss(line1);
        string line2;
        while (getline(ss, line2, '\t'))
        {
            string line3 = istringstream(line2).str();
            rawEntries_line.push_back(line3);
        }
        rawEntries.push_back(rawEntries_line);
    }
    
    // close file stream
    globals.data_fileStream.close();
    
    // if no data then abort
    if (rawEntries.size()==0) {
        cerrAndLog("\nError: data file appears to be empty\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // check same number of values in each row
    int cols = int(rawEntries[0].size());
    for (int i=0; i<int(rawEntries.size()); i++) {
        int row = globals.headerRow_on ? i+2 : i+1;
        if (rawEntries[i].size()!=cols) {
            cerrAndLog("\nError: row "+to_string((long long)row)+string(" in data file does not contain ")+to_string((long long)cols)+string(" entries\n"), globals.outputLog_on, globals.outputLog_fileStream);
            exit(1);
        }
    }
    
    // check that number of columns is compatible with data format options
    int minCols = 2;
    minCols = globals.popCol_on ? minCols+1 : minCols;
    minCols = globals.ploidyCol_on ? minCols+1 : minCols;
    if (cols<minCols) {
        cerrAndLog("\nError: data file must contain at least "+to_string((long long)minCols)+string(" columns under current settings\n"), globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // find which columns in data file correspond to special values
    int pop_startRead = -1;
    int ploidy_startRead = -1;
    int data_startRead = 1;
    if (globals.ploidyCol_on) {
        ploidy_startRead = 1;
        data_startRead++;
    }
    if (globals.popCol_on) {
        pop_startRead = 1;
        if (ploidy_startRead!=-1)
            ploidy_startRead++;
        data_startRead++;
    }
    
    // print these to screen and file
    if (globals.headerRow_on)
        coutAndLog("  row 1 = header line\n", globals.outputLog_on, globals.outputLog_fileStream);
    coutAndLog("  column 1 = individual labels\n", globals.outputLog_on, globals.outputLog_fileStream);
    if (pop_startRead!=-1)
        coutAndLog("  column "+to_string((long long)pop_startRead+1)+" = population of origin\n", globals.outputLog_on, globals.outputLog_fileStream);
    if (ploidy_startRead!=-1)
        coutAndLog("  column "+to_string((long long)ploidy_startRead+1)+" = ploidy\n", globals.outputLog_on, globals.outputLog_fileStream);
    
    // split rawData into special columns and data matrix
    vector<string> indLabels_raw, pop_raw, ploidy_raw;
    vector< vector<string> > rawData;
    vector<string> rawData_line;
    for (int i=0; i<int(rawEntries.size()); i++) {
        rawData_line.clear();
        for (int j=0; j<int(rawEntries[i].size()); j++) {
            
            // read in individual labels
            if (j==0)
                indLabels_raw.push_back(rawEntries[i][j]);
            
            // read in population labels
            if (j==pop_startRead)
                pop_raw.push_back(rawEntries[i][j]);
            
            // read in ploidy for each individual
            if (j==ploidy_startRead)
                ploidy_raw.push_back(rawEntries[i][j]);
            
            // read in data
            if (j>=data_startRead)
                rawData_line.push_back(rawEntries[i][j]);
            
        }
        rawData.push_back(rawData_line);
    }
    
    // finally define number of loci
    globals.loci = int(rawData[0].size());
    
    // convert ploidy labels into integer format
    vector<int> ploidy_raw_int(rawEntries.size(),globals.ploidy);
    if (globals.ploidyCol_on) {
        for (int i=0; i<int(ploidy_raw.size()); i++) {
            double param_double=0;
            istringstream(ploidy_raw[i]) >> param_double;
            if (round(double(param_double))==param_double && round(double(param_double))>=0) {
                ploidy_raw_int[i] = int(param_double);
            } else {
                cerrAndLog("\nError: all ploidy values in data file must be positive integers\n", globals.outputLog_on, globals.outputLog_fileStream);
                exit(1);
            }
        }
    }
    
    // check that ploidy values make sense. Ploidy must be either defined on the first row of an individual (ploidy_format==1), or defined equally for all rows of an individual (ploidy_format==2).
    int thisPloidy=0, countDown=0, ploidy_format=0;
    bool ploidyError = false;
    for (int i=0; i<int(rawEntries.size()); i++) {
        // if new individual
        if (countDown==0) {
            if (ploidy_raw_int[i]!=0) {
                thisPloidy = ploidy_raw_int[i];
                countDown = thisPloidy-1;
            } else {
                ploidyError = true;
            }
        } else {
            if (ploidy_raw_int[i]==0) {
                if (ploidy_format==0 || ploidy_format==1) {
                    ploidy_format = 1;
                    countDown--;
                } else {
                    ploidyError = true;
                }
            } else if (ploidy_raw_int[i]==thisPloidy) {
                if (ploidy_format==0 || ploidy_format==2) {
                    ploidy_format = 2;
                    countDown--;
                } else {
                    ploidyError = true;
                }
            } else {
                ploidyError = true;
            }
        }
        if (ploidyError) {
            cerrAndLog("\nError: when reading ploidy for each individual from file (using the ploidyCol_on=true option) the data file must contain a ploidy value on the first row of each individual, or alternatively the same ploidy value for all rows of an individual\n", globals.outputLog_on, globals.outputLog_fileStream);
            exit(1);
        }
    }
    if (countDown!=0) {
        cerrAndLog("\nError: data file does not contain sufficient rows to account for ploidy values (also check that headerRow_on is specified correctly)\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // similarly for pop values. Pop values must be either defined on the first row of an individual (pop_format==1), or defined equally for all rows of an individual (pop_format==2).
    if (globals.popCol_on) {
        countDown=0;
        string thisPop;
        int pop_format=0;
        bool popError = false;
        for (int i=0; i<int(rawEntries.size()); i++) {
            // if new individual
            if (countDown==0) {
                thisPop = pop_raw[i];
                countDown = ploidy_raw_int[i]-1;
            // if not new individual
            } else {
                if (pop_raw[i]=="") {
                    if (pop_format==0 || pop_format==1) {
                        pop_format = 1;
                        countDown--;
                    } else {
                        popError = true;
                    }
                } else if (pop_raw[i]==thisPop) {
                    if (pop_format==0 || pop_format==2) {
                        pop_format = 2;
                        countDown--;
                    } else {
                        popError = true;
                    }
                } else {
                    popError = true;
                }
            }
            if (popError) {
                cerrAndLog("\nError: when reading population of origin for each individual from file (using the popCol_on=true option) the data file must contain a population value on the first row of each individual, or alternatively the same population value for all rows of an individual\n", globals.outputLog_on, globals.outputLog_fileStream);
                exit(1);
            }
        }
    }
    
    // find all unique alleles, and count number of missing values
    globals.J = vector<int>(globals.loci);
    globals.uniqueAlleles = vector< vector<string> >(globals.loci,vector<string>());
    int missingDataCount = 0;
    for (int j=0; j<globals.loci; j++) {
        for (int i=0; i<int(rawData.size()); i++) {
            if (rawData[i][j]==globals.missingData) {
                missingDataCount++;
            } else {
                if (find(globals.uniqueAlleles[j].begin(), globals.uniqueAlleles[j].end(), rawData[i][j]) == globals.uniqueAlleles[j].end()) {
                    globals.uniqueAlleles[j].push_back(rawData[i][j]);
                }
            }
        }
        globals.J[j] = int(globals.uniqueAlleles[j].size());
    }
    
    // reformat data into simple list: data[individual][locus][gene copy]. Recode values to simple list of integers, with 0 meaning missing data. Keep record of missing data elements in data_missing (1=missing, 0=present). Store indLabels_vec, pop_vec, ploidy_vec and missing_vec values.
    countDown=0;
    int ind=-1, gene_copy=0;
    for (int i=0; i<int(rawEntries.size()); i++) {
        // if new individual
        if (countDown==0) {
            ind++;
            globals.indLabels_vec.push_back(indLabels_raw[i]);
            if (globals.popCol_on) {
                globals.pop_vec.push_back(pop_raw[i]);
            } else {
                globals.pop_vec.push_back("NA");
            }
            globals.ploidy_vec.push_back(ploidy_raw_int[i]);
            globals.missing_vec.push_back(0);
            countDown = globals.ploidy_vec[ind]-1;
            gene_copy = 0;
            
            globals.data.push_back(vector< vector<int> >(globals.loci,vector<int>(globals.ploidy_vec[ind])));
            for (int l=0; l<globals.loci; l++) {
                if (rawData[i][l]==globals.missingData) {
                    globals.data[ind][l][gene_copy] = 0;
                    globals.missing_vec[ind]++;
                } else {
                    globals.data[ind][l][gene_copy] = whichFirst(globals.uniqueAlleles[l],rawData[i][l])+1;
                }
            }
        // if not new individual
        } else {
            countDown--;
            gene_copy++;
            for (int l=0; l<globals.loci; l++) {
                if (rawData[i][l]==globals.missingData) {
                    globals.data[ind][l][gene_copy] = 0;
                    globals.missing_vec[ind]++;
                } else {
                    globals.data[ind][l][gene_copy] = whichFirst(globals.uniqueAlleles[l],rawData[i][l])+1;
                }
            }
        }
    }
    
    // finally define number of individuals and gene copies
    globals.n = int(globals.indLabels_vec.size());
    globals.geneCopies = sum(globals.ploidy_vec)*globals.loci;
    
    // check that there is at least one allele at every locus
    for (int l=0; l<globals.loci; l++) {
        bool tmp = true;
        for (int i=0; i<globals.n; i++) {
            for (int p=0; p<globals.ploidy_vec[i]; p++) {
                if (globals.data[i][l][p]!=0)
                    tmp = false;
            }
        }
        if (tmp) {
            coutAndLog("\nError: there must be at least one non-missing observation at every locus\n", globals.outputLog_on, globals.outputLog_fileStream);
            exit(1);
        }
    }
    
    // find number of unique populations and store pop index of each individual
    if (globals.popCol_on) {
        globals.pop_index = vector<int>(globals.n);
        for (int i=0; i<globals.n; i++) {
            if (find(globals.uniquePops.begin(), globals.uniquePops.end(), globals.pop_vec[i]) == globals.uniquePops.end()) {
                globals.uniquePops.push_back(globals.pop_vec[i]);
                globals.pop_index[i] = int(globals.uniquePops.size())-1;
            } else {
                globals.pop_index[i] = int(find(globals.uniquePops.begin(), globals.uniquePops.end(), globals.pop_vec[i]) - globals.uniquePops.begin());
            }
        }
        globals.uniquePop_counts = vector<int>(globals.uniquePops.size());
        for (int i=0; i<globals.n; i++) {
            globals.uniquePop_counts[globals.pop_index[i]]++;
        }
    }
    
    // print basic properties to screen and file
    if (globals.popCol_on)
        coutAndLog("  unique populations = "+to_string((long long)globals.uniquePops.size())+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
    
    coutAndLog("  individuals = "+to_string((long long)globals.n)+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
    coutAndLog("  loci = "+to_string((long long)globals.loci)+string("\n"), globals.outputLog_on, globals.outputLog_fileStream);
    
    stringstream alleles;
    alleles << "  alleles per locus = {" << globals.J[0];
    for (int l=1; l<globals.loci; l++) {
        alleles << "," << globals.J[l];
    }
    alleles << "}\n";
    coutAndLog(alleles.str(), globals.outputLog_on, globals.outputLog_fileStream);
    
    coutAndLog("  missing observations = "+to_string((long long)missingDataCount)+string(" of ")+to_string((long long)sum(globals.ploidy_vec)*globals.loci)+string("\n\n"), globals.outputLog_on, globals.outputLog_fileStream);    
    
}

//------------------------------------------------
// check that the options chosen (in terms of data and parameters) make sense
void checkOptions(globals &globals) {
    
    /*
     The following checks are carried out (in this order):
     - fixLabels_on=true if producing any Qmatrix output
     - popCol_on=true if producing Qmatrix_pop files
     - admix_on=true if
        - producing Qmatrix_gene files
        - producing ML admixture freqs
     - outputQmatrix_?_on=true if outputQmatrixError_?_on=true
     - mainRepeats>1 if producing QmatrixError files
     - EMalgorithm_on=true if
        - producing model comp stats
        - producing ML output
     - (Kmax-Kmin)>=2 and mainRepeats>1 if outputEvanno_on=true
     - if using exhaustive approach, check that number of partitions required to sum over is reasonable. Otherwise throw warning
     - ensure that alpha and alphaPropSD vectors are of length 1 or (Kmax-Kmin+1). If former then duplicate value as needed.
     */
    
    // force fixLabels_on=true if producing any Qmatrix output
    if (!globals.fixLabels_on && (globals.outputQmatrix_ind_on || globals.outputQmatrix_pop_on || globals.outputQmatrix_gene_on)) {
        cerrAndLog("\nError: the label switching problem must be solved in order to produce Qmatrix output. Either set fixLabels_on to true, or turn off all Qmatrix output.\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // force popCol_on=true if producing Qmatrix_pop files
    if (!globals.popCol_on && globals.outputQmatrix_pop_on) {
        cerrAndLog("\nError: unable to produce population level Qmatrix, as population data not currently being used. Either read in population data (by setting popCol_on to true), or stop producing population level Qmatrix output (by setting outputQmatrix_pop_on to false and outputQmatrixError_pop_on to false).\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // force admix_on=true if producing Qmatrix_gene files
    if (!globals.admix_on && globals.outputQmatrix_gene_on) {
        cerrAndLog("\nError: admixture model must be turned on in order to produce gene level Qmatrix. Either switch to the admixture model (by setting admix_on to true) or stop producing gene level Qmatrix output (by setting outputQmatrix_gene_on and outputQmatrixError_gene_on to false).\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // force admix_on=true if producing maximum likelihood admixture frequencies
    if (!globals.admix_on && globals.outputMaxLike_admixFreqs_on) {
        cerrAndLog("\nError: admixture model must be turned on in order to produce maximum likelihood admixture frequency output. Either switch to the admixture model (by setting admix_on to true) or stop producing maximum likelihood admixture frequency output (by setting outputMaxLike_admixFreqs_on to false).\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // force outputQmatrix_?_on=true if outputQmatrixError_?_on=true for all Qmatrix types
    if ((globals.outputQmatrixError_ind_on && !globals.outputQmatrix_ind_on) || (globals.outputQmatrixError_pop_on && !globals.outputQmatrix_pop_on) || (globals.outputQmatrixError_gene_on && !globals.outputQmatrix_gene_on)) {
        cerrAndLog("\nError: if outputQmatrixError options are turned on then the corresponding outputQmatrix options must also be turned on.\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // force mainRepeats>1 if outputQmatrixError_?_on=true for all Qmatrix types
    if (globals.mainRepeats==1 && (globals.outputQmatrixError_ind_on || globals.outputQmatrixError_pop_on || globals.outputQmatrixError_gene_on)) {
        cerrAndLog("\nError: mainRepeats must be greater than 1 when any outputQmatrixError options are turned on.\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // force EMalgorithm_on=true if producing model comparison statistics
    if (!globals.EMalgorithm_on && globals.outputComparisonStatistics_on) {
        cerrAndLog("\nError: EM algorithm must be turned on in order to compute model comparison statistics. Either turn on the EM algorithm (by setting EMalgorithm_on to true) or stop producing model comparison statistics (by setting outputComparisonStatistics_on to false).\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // force EMalgorithm_on=true if producing maximum likelihood output
    if (!globals.EMalgorithm_on && (globals.outputMaxLike_alleleFreqs_on || globals.outputMaxLike_admixFreqs_on)) {
        cerrAndLog("\nError: EM algorithm must be turned on in order to produce maximum likelihood output files. Either turn on the EM algorithm (by setting EMalgorithm_on to true) or stop producing maximum likelihood output (by setting outputMaxLike_alleleFreqs_on and outputMaxLike_admixFreqs_on to false).\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // force (Kmax-Kmin)>=2 and mainRepeats>1 if outputEvanno_on=true
    if ((globals.Kmax-globals.Kmin)<2 && globals.outputEvanno_on) {
        cerrAndLog("\nError: Kmax must be at least 2 greater than Kmin when calculating Evanno's delta K.\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    if (globals.mainRepeats==1 && globals.outputEvanno_on) {
        cerrAndLog("\nError: mainRepeats must be greater than 1 when calculating Evanno's delta K.\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    // if using exhaustive approach, check that number of partitions required to sum over is reasonable. Otherwise throw warning
    if (globals.exhaustive_on && !globals.suppressWarning1_on) {
        
        // calculate number of partitions. Number of ways of splitting n elements into i groups is given by the Stirling number of the second kind {n, i}. Number of ways of splitting n elements into UP TO K groups is given by sum_{i=1}^k {n, i}. Calculate this for all K that will be explored under current settings. The total number of partitions (in log space) is given by log_partitions.
        double log_partitions = log(double(0));
        for (int K=globals.Kmin; K<=globals.Kmax; K++) {
            vector<double> StirlingNumbers = globals.admix_on ? StirlingSecond(sum(globals.ploidy_vec)*globals.loci) : StirlingSecond(globals.n);
            for (int i=0; (i<K) & (i<int(StirlingNumbers.size())); i++) {
                log_partitions = logSum(log_partitions,StirlingNumbers[i]);
            }
        }
        
        // if very large then error
        if (log_partitions>180) {
            cerrAndLog("\nWarning: for this combination of data and parameters the number of partitions explored in the exhaustive approach is greater than the number of protons in the observable universe, and the program is unlikely to finish within your lifetime. Either turn off the exhaustive approach, or to override this warning and continue with the analysis insert the optional parameter \"suppressWarning1_on\" with value \"true\" in the parameters file.\n", globals.outputLog_on, globals.outputLog_fileStream);
            exit(1);
        }
        if (log_partitions>log(double(1e7))) {
            cerrAndLog("\nWarning: for this combination of data and parameters the number of partitions explored in the exhaustive approach is extremely large (greater than 1e7 in total), and the program may take a very long time to execute. Either turn off the exhaustive approach, or to override this warning and continue with the analysis insert the optional parameter \"suppressWarning1_on\" with value \"true\" in the parameters file.\n", globals.outputLog_on, globals.outputLog_fileStream);
            exit(1);
        }
        
    }
    
    // ensure that alpha and alphaPropSD vectors are of length 1 or (Kmax-Kmin+1). If former then duplicate value as needed.
    if (int(globals.alpha.size())==1) {
        globals.alpha = vector<double>(globals.Kmax-globals.Kmin+1, globals.alpha[0]);
    } else if (int(globals.alpha.size())!=(globals.Kmax-globals.Kmin+1)) {
        cerrAndLog("\nWarning: alpha must contain either a single value to apply to all K, or a comma-separated list of values of length (Kmax-Kmin+1) to apply to each K individually.\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
    if (int(globals.alphaPropSD.size())==1) {
        globals.alphaPropSD = vector<double>(globals.Kmax-globals.Kmin+1, globals.alphaPropSD[0]);
    } else if (int(globals.alphaPropSD.size())!=(globals.Kmax-globals.Kmin+1)) {
        cerrAndLog("\nWarning: alphaPropSD must contain either a single value to apply to all K, or a comma-separated list of values of length (Kmax-Kmin+1).\n", globals.outputLog_on, globals.outputLog_fileStream);
        exit(1);
    }
    
}
