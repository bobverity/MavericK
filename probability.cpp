//
//  MavericK
//  probability.cpp
//
//  Created: Bob on 22/09/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "probability.h"

using namespace std;

//-- set random seed --
random_device rd;
default_random_engine generator(rd());
uniform_real_distribution<double> uniform_0_1(0.0,1.0);

//------------------------------------------------
// draw from uniform(a,b) distribution
double runif1(double a, double b) {
    uniform_real_distribution<double> uniform_a_b(a,b);
    return(uniform_a_b(generator));
}

//------------------------------------------------
// draw from gamma(shape,rate) distribution
double rgamma1(double shape, double rate) {
    gamma_distribution<double> rgamma(shape,1/rate);
	double x = rgamma(generator);

	// check for zero or infinite values (catches bug present in Visual Studio 2010)
	if (x==0)
		x = UNDERFLO;
	while ((1.0/x)==0)
		x = rgamma(generator);

    return(x);
}

//------------------------------------------------
// sample from given probability vector that sums to pSum
int sample1(vector<double> &p, double pSum) {
    double rand = pSum*uniform_0_1(generator);
    double z = 0;
    for (int i=0; i<int(p.size()); i++) {
        z += p[i];
        if (rand<z)
            return i+1;
    }
    return(0);
}

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd) {
    normal_distribution<double> rnorm(mean,sd);
    return(rnorm(generator));
}
