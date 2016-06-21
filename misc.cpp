//
//  MavericK
//  misc.cpp
//
//  Created: Bob on 22/09/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "misc.h"

using namespace std;

//------------------------------------------------
// define very small number for catching underflow problems
// DEFINED IN HEADER

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double logSum(double logA, double logB) {
    if (logA-logB > 100) {
        return(logA);
    } else if (logB-logA > 100) {
        return(logB);
    }
    double output = (logA<logB) ? logB + log(1+exp(logA-logB)) : logA + log(1+exp(logB-logA));
    return(output);
}

//------------------------------------------------
// mean of vector (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// sample or population variance of vector (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// calculate harmonic mean of vector, accounting for underflow. x must be in log space.
double harmonicMean(vector<double> &x) {
    double output = 0;
    double minVal =  *min_element(begin(x),end(x));
    for (int i=0; i<int(x.size()); i++)
        output += 1.0/exp(x[i]-minVal);
    output = minVal- log(output/x.size());
    return(output);
}

//------------------------------------------------
// exponentiate and normalise vector to sum to 1
vector<double> normalise_log(vector<double> x) {
    
    vector<double> output(x.size());
    double maxVal = *max_element(begin(x),end(x));
    double sumVal = 0;
    for (int i=0; i<int(x.size()); i++) {
        output[i] = exp(x[i]-maxVal);
        sumVal += output[i];
    }
    for (int i=0; i<int(x.size()); i++) {
        output[i] /= sumVal;
    }
    return(output);
}

//------------------------------------------------
// exponentiate and normalise random variables to sum to 1 by simulation. Return mean values, and upper and lower 95% confidence intervals (q0.025 and q0.975).
void normalise_log_sim(vector<double> &normMean, vector<double> &LL, vector<double> &UL, vector<double> x_mean, vector<double> x_sd, int draws) {
    
    // subtract maximum value from x_mean (this has no effect on final outcome but prevents under/overflow)
    int n = int(x_mean.size());
    double maxVal = *max_element(begin(x_mean),end(x_mean));
    for (int j=0; j<n; j++) {
        x_mean[j] -= maxVal;
    }
    
    // draw random values of Z (exponentiated and normalised draws)
    vector<double> X(n);
    vector<double> Y(n);
    double Ysum;
    vector< vector<double> > Z(n,vector<double>(draws));
    for (int i=0; i<draws; i++) {
        Ysum = 0;
        for (int j=0; j<n; j++) {
            X[j] = rnorm1(x_mean[j], x_sd[j]);
            Y[j] = exp(X[j]);
            Ysum += Y[j];
        }
        for (int j=0; j<n; j++) {
            Z[j][i] = Y[j]/Ysum;
            if (Ysum==0)
                Z[j][i] = 1/double(n);
        }
    }
    
    // sort values of Z
    for (int j=0; j<n; j++) {
        sort(Z[j].begin(),Z[j].end());
    }
    
    // store mean and 95% CI values
    for (int j=0; j<int(X.size()); j++) {
        //median[j] = Z[j][round(draws*0.5)];
        normMean[j] = mean(Z[j]);
        LL[j] = Z[j][round(draws*0.025)];
        UL[j] = Z[j][round(draws*0.975)];
    }
    
}

//------------------------------------------------
// return unique elements in a vector (templated for different data types)
// DEFINED IN HEADER


//------------------------------------------------
// return position of FIRST element in vector 'haystack' equal to 'needle' (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// return vector of first n Stirling numbers of the first kind
vector<int> StirlingFirst(int n) {
    
    vector<int> oldVec(2,1);
    vector<int> newVec(2);
    
    // simple answer if n<=0, n==1 or n==2
    if (n<=0)
        return(vector<int>(1));
    else if (n==1)
        return(vector<int>(1,1));
    else if (n==2)
        return(oldVec);
    
    // otherwise calculate by looping
    for (int i=2; i<n; i++) {
        newVec[0] = oldVec[0]*i;
        for (int j=1; j<int(oldVec.size()); j++) {
            newVec[j] = oldVec[j]*i + oldVec[j-1];
        }
        newVec.push_back(1);
        oldVec = vector<int>(newVec);
    }
    
    return(newVec);
}

//------------------------------------------------
// return vector of first n Stirling numbers of the second kind (calculate and output in log space)
vector<double> StirlingSecond(int n) {
    
    vector<double> oldVec(2,0);
    vector<double> newVec(2,0);
    
    // simple answer if n<=0, n==1 or n==2
    if (n<=0)
        return(vector<double>(1,log(0.0)));
    else if (n==1)
        return(vector<double>(1));
    else if (n==2)
        return(oldVec);
    
    // otherwise calculate by looping
    for (int i=2; i<n; i++) {
        for (int j=1; j<int(oldVec.size()); j++) {
            newVec[j] = logSum(oldVec[j-1],oldVec[j]+log(double(j+1)));
        }
        newVec.push_back(0);
        oldVec = vector<double>(newVec);
    }
    
    return(newVec);
}

//------------------------------------------------
// exit when user presses enter
void pauseExit() {
    cin.sync();
    cin.get();
    exit(1);
}

//------------------------------------------------
// helper function for printing a single value (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// same as printVector, but exponentiates values prior to printing (overloaded for int and double)
void printVector_exp(vector<int> &x) {
    for (int i=0; i<int(x.size()); i++) {
        cout << exp(double(x[i])) << " ";
    }
    cout << "\n";
}
void printVector_exp(vector<double> &x) {
    for (int i=0; i<int(x.size()); i++) {
        cout << exp(x[i]) << " ";
    }
    cout << "\n";
}

//------------------------------------------------
// same as printMatrix, but exponentiates values prior to printing (overloaded for int and double)
void printMatrix_exp(vector< vector<int> > &M) {
    for (int i=0; i<int(M.size()); i++) {
        for (int j=0; j<int(M[i].size()); j++) {
            cout << exp(double(M[i][j])) << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}
void printMatrix_exp(vector< vector<double> > &M) {
    for (int i=0; i<int(M.size()); i++) {
        for (int j=0; j<int(M[i].size()); j++) {
            cout << exp(M[i][j]) << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

//------------------------------------------------
// same as printArray, but exponentiates values prior to printing (overloaded for int and double)
void printArray_exp(vector< vector< vector<int> > > &x) {
    for (int i=0; i<int(x.size()); i++) {
        cout << "--- slice " << i+1 << " ---\n";
        for (int j=0; j<int(x[i].size()); j++) {
            for (int k=0; k<int(x[i][j].size()); k++) {
                cout << exp(double(x[i][j][k])) << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }
    cout << "\n";
}
void printArray_exp(vector< vector< vector<double> > > &x) {
    for (int i=0; i<int(x.size()); i++) {
        cout << "--- slice " << i+1 << " ---\n";
        for (int j=0; j<int(x[i].size()); j++) {
            for (int k=0; k<int(x[i][j].size()); k++) {
                cout << exp(x[i][j][k]) << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }
    cout << "\n";
}

//------------------------------------------------
// check that the value in param_string is a Boolean value, and drop result into param_final. If fail test then print error message to screen, and to logFileStream if writeErrorToFile is true. param_name gives the name of the parameter to be used in error message.
void checkBoolean(string param_string, bool &param_final, string param_name, bool writeErrorToFile, ofstream &logFileStream) {
    if (param_string=="1" || param_string=="1\r" || param_string=="true" || param_string=="true\r" || param_string=="True" || param_string=="True\r" ||  param_string=="TRUE" || param_string=="TRUE\r" || param_string=="t" || param_string=="t\r" || param_string=="T" || param_string=="T\r") {
        param_final = true;
    } else if (param_string=="0" || param_string=="0\r" || param_string=="false" || param_string=="false\r" || param_string=="False" || param_string=="False\r" || param_string=="FALSE" || param_string=="FALSE\r" || param_string=="f" || param_string=="f\r" || param_string=="F" || param_string=="F\r") {
        param_final = false;
    } else {
        cerrAndLog("\nError: '"+param_name+"' parameter must take on a Boolean value (either 1/t/true/T/True/TRUE or 0/f/false/F/False/FALSE)\n", writeErrorToFile, logFileStream);
        exit(1);
    }
}

//------------------------------------------------
// check that the value in param_string is an integer value, and drop result into param_final. If fail test then print error message to screen, and to logFileStream if writeErrorToFile is true. param_name gives the name of the parameter to be used in error message.
void checkInteger(string param_string, int &param_final, string param_name, bool writeErrorToFile, ofstream &logFileStream) {
    double param_double;
    istringstream(param_string) >> param_double;
    if (round(double(param_double))==param_double) {
        param_final = int(param_double);
    } else {
        cerrAndLog("\nError: '"+param_name+"' parameter must be an integer\n", writeErrorToFile, logFileStream);
        exit(1);
    }
}

//------------------------------------------------
// check that param_value is greater than zero (overloaded for integer and double). If fail test then print error message to screen, and to logFileStream if writeErrorToFile is true. param_name gives the name of the parameter to be used in error message.
void checkGrZero(string param_name, int param_value, bool writeErrorToFile, ofstream &logFileStream) {
    if (!(param_value>0)) {
        cerrAndLog("\nError: '"+param_name+"' parameter must be greater than zero\n", writeErrorToFile, logFileStream);
        exit(1);
    }
}
void checkGrZero(string param_name, double param_value, bool writeErrorToFile, ofstream &logFileStream) {
    if (!(param_value>0)) {
        cerrAndLog("\nError: '"+param_name+"' parameter must be greater than zero\n", writeErrorToFile, logFileStream);
        exit(1);
    }
}

//------------------------------------------------
// check that param_value is greater than or equal to zero (overloaded for integer and double). If fail test then print error message to screen, and to logFileStream if writeErrorToFile is true. param_name gives the name of the parameter to be used in error message.
void checkGrEqZero(string param_name, int param_value, bool writeErrorToFile, ofstream &logFileStream) {
    if (!(param_value>=0)) {
        cerrAndLog("\nError: '"+param_name+"' parameter must be greater than or equal to zero\n", writeErrorToFile, logFileStream);
        exit(1);
    }
}
void checkGrEqZero(string param_name, double param_value, bool writeErrorToFile, ofstream &logFileStream) {
    if (!(param_value>=0)) {
        cerrAndLog("\nError: '"+param_name+"' parameter must be greater than or equal to zero\n", writeErrorToFile, logFileStream);
        exit(1);
    }
}

//------------------------------------------------
// increment grouping according to restricted growth constraint. For example, {1,1,2,1,2,3} increments to {1,1,2,1,3,1}, making it possible to list all unique partitions. Maximum value is limited to maxVal, for example {1,1,2,1,2,2} increments to {1,1,2,2,1,1} when maxVal=2. Returns group[0]=0 when no more partitions.
void increment_restrictedGrowth(vector<int> &group, int maxVal) {
    
    // special escape condition if maxVal==1
    if (maxVal==1) {
        group[0] = 0;
        return;
    }
    
    // starting at the far right end, keep incrementing individual elements until restricted growth function condition is met.
    int target = int(group.size());
    group[target-1] ++;
    int targetMax = *max_element(begin(group),begin(group)+target-1);
    while ((group[target-1]>(targetMax+1)) || group[target-1]>maxVal) {
        group[target-1] = 1;
        target --;
        group[target-1] ++;
        targetMax = *max_element(begin(group),begin(group)+target-1);
    }
    if (target==1) {
        group[0] = 0;
    }
}

//------------------------------------------------
// return double value as string, replacing nan values with "NA"
string process_nan(double x) {
    string s = to_string((long double)x);
    if (isnan(x))
        s = "NA";
    return(s);
}

//------------------------------------------------
// calculate total autocorrelation on a vector of values, where total autocorrelation is defined as 1+2*(sum over lags until reach zero). This number is equivalent to the number of iterations needed to obtain one approximately independent draw (1 being the best). Note that this function calculates raw autocorrelation of the vector in increasingly large jumps, leading to fine resolution near lag1 and coarser resolution as we get further away.
double calculateAutoCorr(vector<double> &v) {
    
    // get basic properties of vector
    int v_size = int(v.size());
    double mu = mean(v);
    double sigma2 = var(v);
    
    // calculate x range over which autocorrelation will be calculated
    vector<int> x;
    int x_largest = 0;
    int dx = 1;
    bool break_on = false;
    for (int l=0; l<100; l++) {
        for (int j=0; j<100; j++) {
            x_largest += dx;
            if (x_largest>v_size) {
                break_on = true;
                break;
            }
            x.push_back(x_largest);
        }
        if (break_on) {
            break;
        }
        dx *= 2;
    }
    
    // calculate raw autocorrelation
    int x_size = int(x.size());
    vector<double> yy(x_size);
    vector<double> n(x_size);
    for (int i=0; i<v_size; i++) {
        for (int j=0; j<x_size; j++) {
            if ((i+x[j])>(v_size-1)) {
                continue;
            }
            yy[j] += (v[i]-mu)*(v[i+x[j]]-mu)/sigma2;
            n[j]++;
        }
    }
    for (int j=0; j<x_size; j++) {
        yy[j] /= double(n[j]);
    }
    
    // calculate total autocorrelation
    double autoCorr = yy[0];
    for (int i=1; i<x_size; i++) {
        dx = x[i]-x[i-1];
        if (yy[i]>0) {
            autoCorr += dx*yy[i-1] - 0.5*(dx+1)*(yy[i-1]-yy[i]);
        } else {
            dx = floor((dx*yy[i-1])/(yy[i-1]-yy[i]));
            autoCorr += dx*yy[i-1] - 0.5*(dx+1)*(yy[i-1]-yy[i]);
            break;
        }
    }
    autoCorr = 1 + 2*autoCorr;
    
    // correct for possible issues
    if (sigma2==0 || autoCorr<1)
        autoCorr = 1;
    
    return(autoCorr);
}
