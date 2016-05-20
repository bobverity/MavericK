//
//  MavericK
//  misc.h
//
//  Created: Bob on 22/09/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Miscellaneous functions.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__misc__
#define __Maverick1_0__misc__

#include <iostream>
#include <random>
#include <algorithm>

#include "probability.h"
#include "readIn.h"
#include "OSfunctions.h"

//------------------------------------------------
// define very small number for catching underflow problems
#define UNDERFLO   1e-100

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types).
template<class TYPE>
TYPE sum(std::vector<TYPE> &x) {
    TYPE output = 0;
    for (int i=0; i<int(x.size()); i++)
        output += x[i];
    return(output);
}

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double logSum(double logA, double logB);

//------------------------------------------------
// mean of vector (templated for different data types)
template<class TYPE>
double mean(std::vector<TYPE> &x) {
    return(sum(x)/double(x.size()));
}

//------------------------------------------------
// sample or population variance of vector (templated for different data types)
template<class TYPE>
double var(std::vector<TYPE> &x, double mu=0, bool sampleMean=true) {
    int n = int(x.size());
    if (n==1)
        return(0);
    double xSum=0, xSumSquare=0;
    for (int i=0; i<n; i++) {
        xSum += x[i];
        xSumSquare += x[i]*x[i];
    }
    double output = sampleMean ? (xSumSquare-xSum*xSum/n)/(n-1) : (xSumSquare-2*mu*xSum+n*mu*mu)/n;
    output = (output<0) ? 0 : output;
    return(output);
}

//------------------------------------------------
// calculate harmonic mean of vector, accounting for underflow. x must be in log space.
double harmonicMean(std::vector<double> &x);

//------------------------------------------------
// exponentiate and normalise vector to sum to 1
std::vector<double> normalise_log(std::vector<double> x);

//------------------------------------------------
// exponentiate and normalise random variables to sum to 1 by simulation. Return median values, and upper and lower 95% confidence intervals (q0.025 and q0.975).
void normalise_log_sim(std::vector<double> &median, std::vector<double> &LL, std::vector<double> &UL, std::vector<double> x_mean, std::vector<double> x_sd, int draws);

//------------------------------------------------
// return unique elements in a vector (templated for different data types)
template<class TYPE>
std::vector<TYPE> uniques(std::vector<TYPE> x) {
    std::vector<TYPE> output;
    for (int i=0; i<x.size(); i++) {
        if (find(output.begin(), output.end(), x[i]) == output.end()) {
            output.push_back(x[i]);
        }
    }
    return(output);
}

//------------------------------------------------
// return position of FIRST element in vector 'haystack' equal to 'needle' (templated for different data types)
template<class TYPE>
int whichFirst(std::vector<TYPE> haystack, TYPE needle) {
    int index = 0;
    for (int i=0; i<int(haystack.size()); i++) {
        if (haystack[i]==needle) {
            break;
        }
        index++;
    }
    return(index);
}

//------------------------------------------------
// return vector of first n Stirling numbers of the first kind
std::vector<int> StirlingFirst(int n);

//------------------------------------------------
// return vector of first n Stirling numbers of the second kind (calculate and output in log space)
std::vector<double> StirlingSecond(int n);

//------------------------------------------------
// exit when user presses enter
void pauseExit();

//------------------------------------------------
// helper function for printing a single value (templated for different data types)
template<class TYPE>
void print(TYPE x) {
    std::cout << x << "\n";
}

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
template<class TYPE>
void printVector(std::vector<TYPE> &x) {
    for (int i=0; i<x.size(); i++) {
        std::cout << x[i] << " ";
    }
    std::cout << "\n";
}

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
template<class TYPE>
void printMatrix(std::vector< std::vector<TYPE> > &M) {
    for (int i=0; i<M.size(); i++) {
        for (int j=0; j<M[i].size(); j++) {
            std::cout << M[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
template<class TYPE>
void printArray(std::vector< std::vector< std::vector<TYPE> > > &x) {
    for (int i=0; i<x.size(); i++) {
        std::cout << "--- slice " << i+1 << " ---\n";
        for (int j=0; j<x[i].size(); j++) {
            for (int k=0; k<x[i][j].size(); k++) {
                std::cout << x[i][j][k] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

//------------------------------------------------
// same as printVector, but exponentiates values prior to printing (overloaded for int and double)
void printVector_exp(std::vector<int> &x);
void printVector_exp(std::vector<double> &x);

//------------------------------------------------
// same as printMatrix, but exponentiates values prior to printing (overloaded for int and double)
void printMatrix_exp(std::vector< std::vector<int> > &M);
void printMatrix_exp(std::vector< std::vector<double> > &M);

//------------------------------------------------
// same as printArray, but exponentiates values prior to printing (overloaded for int and double)
void printArray_exp(std::vector< std::vector< std::vector<int> > > &x);
void printArray_exp(std::vector< std::vector< std::vector<double> > > &x);

//------------------------------------------------
// check that the value in param_string is a Boolean value, and drop result into param_final. If fail test then print error message to screen, and to logFileStream if writeErrorToFile is true. param_name gives the name of the parameter to be used in error message.
void checkBoolean(std::string param_string, bool &param_final, std::string param_name, bool writeErrorToFile, std::ofstream &logFileStream);

//------------------------------------------------
// check that the value in param_string is an integer value, and drop result into param_final. If fail test then print error message to screen, and to logFileStream if writeErrorToFile is true. param_name gives the name of the parameter to be used in error message.
void checkInteger(std::string param_string, int &param_final, std::string param_name, bool writeErrorToFile, std::ofstream &logFileStream);

//------------------------------------------------
// check that param_value is greater than zero (overloaded for integer and double). If fail test then print error message to screen, and to logFileStream if writeErrorToFile is true. param_name gives the name of the parameter to be used in error message.
void checkGrZero(std::string param_name, int param_value, bool writeErrorToFile, std::ofstream &logFileStream);
void checkGrZero(std::string param_name, double param_value, bool writeErrorToFile, std::ofstream &logFileStream);

//------------------------------------------------
// check that param_value is greater than or equal to zero (overloaded for integer and double). If fail test then print error message to screen, and to logFileStream if writeErrorToFile is true. param_name gives the name of the parameter to be used in error message.
void checkGrEqZero(std::string param_name, int param_value, bool writeErrorToFile, std::ofstream &logFileStream);
void checkGrEqZero(std::string param_name, double param_value, bool writeErrorToFile, std::ofstream &logFileStream);

//------------------------------------------------
// increment grouping according to restricted growth constraint. For example, {1,1,2,1,2,3} increments to {1,1,2,1,3,1}, making it possible to list all unique partitions. Maximum value is limited to maxVal, for example {1,1,2,1,2,2} increments to {1,1,2,2,1,1} when maxVal=2. Returns group[0]=0 when no more partitions.
void increment_restrictedGrowth(std::vector<int> &group, int maxVal);

//------------------------------------------------
// return double value as string, replacing nan values with "NA"
std::string process_nan(double x);

//------------------------------------------------
// calculate total autocorrelation on a vector of values, where total autocorrelation is defined as 1+2*(sum over lags until reach zero). This number is equivalent to the number of iterations needed to obtain one approximately independent draw (1 being the best). Note that this function calculates raw autocorrelation of the vector in increasingly large jumps, leading to fine resolution near lag1 and coarser resolution as we get further away.
double calculateAutoCorr(std::vector<double> &v);

#endif
