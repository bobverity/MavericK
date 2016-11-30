//
//  MavericK
//  readIn.h
//
//  Created: Bob on 24/09/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Functions for reading in data and parameters files, for parsing command line arguments, and for processing the chosen parameter set to make sure that it makes sense.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__readIn__
#define __Maverick1_0__readIn__

#include <iostream>
#include <sstream>
#include "globals.h"
#include "misc.h"
#include "OSfunctions.h"

//------------------------------------------------
// write to file if condition is true
void writeToFile(std::string message, bool writeToFile, std::ofstream &fileStream);

//------------------------------------------------
// write message to screen, and to file if writeToFile is true
void coutAndLog(std::string message, bool writeToFile, std::ofstream &fileStream);

//------------------------------------------------
// write error to screen, and to file if writeToFile is true
void cerrAndLog(std::string message, bool writeToFile, std::ofstream &logFileStream);

//------------------------------------------------
// safely open input file stream, otherwise error. Option to write errors to file if writeErrorToFile is true.
std::ifstream safe_ifstream(std::string fileName, bool writeToFile, std::ofstream &logFileStream);

//------------------------------------------------
// safely open output file stream, otherwise error. Option to write errors to file if writeErrorToFile is true.
std::ofstream safe_ofstream(std::string fileName, bool writeToFile, std::ofstream &logFileStream);

//------------------------------------------------
// convert string to double if possible, otherwise error
double stringToDouble(std::string s, bool writeToFile, std::ofstream &logFileStream);

//------------------------------------------------
// convert double to int (if no remainder), otherwise error
int doubleToInt(double x, bool writeToFile, std::ofstream &logFileStream);

//------------------------------------------------
// equivalent to getline() function, but searches for '\r', '\n' or '\r\n' (therefore Max, PC, Linux compatible)
std::istream& safe_getline(std::istream& is, std::string& t);

//------------------------------------------------
// if argument argv[i] is equal to argName then read item [i+1] into savePath
void readPath(std::string argName, std::string &savePath, int argc, const char * argv[], int i);

//------------------------------------------------
// read in parameters as strings
void readParameters(globals &globals);

//------------------------------------------------
// if argument argv[i] is equal to argName then read item [i+1] into argValue
void readArgument(std::string argName, globals &globals, int argc, const char * argv[], int i);

//------------------------------------------------
// read in command line arguments as strings
void readCommandLine(globals &globals, int argc, const char * argv[]);

//------------------------------------------------
// check parameters, stored in pair<string,int> format. The input i filters the second element in this pair.
void checkParameters(globals &globals, int i);

//------------------------------------------------
// read in data file
void readData(globals &globals);

//------------------------------------------------
// check that the options chosen (in terms of data and parameters) make sense
void checkOptions(globals &globals);

#endif
