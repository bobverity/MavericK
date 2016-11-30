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

extern int glob1;
extern int glob2;

//-- set random seed --
random_device rd;
default_random_engine generator(rd());
//default_random_engine generator(1);

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
// draw from beta(shape1,shape2) distribution
double rbeta1(double shape1, double shape2) {
    double x1 = rgamma1(shape1,1.0);
    double x2 = rgamma1(shape2,1.0);
    return(x1/double(x1+x2));
}

//------------------------------------------------
// sample from given probability vector that sums to pSum
int sample1(vector<double> &p, double pSum) {
    const double rand = pSum*uniform_0_1(generator);
    double z = 0;
    for (unsigned int i=0; i<int(p.size()); i++) {
        z += p[i];
        if (rand<z)
            return i+1;
    }
    return(0);
}

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal probability
int sample2(int a, int b) {
    int z = floor(runif1(a, b+1));
    return(z);
}

//------------------------------------------------
// sample n terms without replacement from vector x
vector<int> sample3(vector<int> x, int n) {
    int rand1, dummy;
    int xSize = int(x.size());
    for (int i=0; i<n; i++) {
        rand1 = sample2(i,xSize-1);
        dummy = x[i];
        x[i] = x[rand1];
        x[rand1] = dummy;
    }
    vector<int> output(x.begin(),x.begin()+n);
    return(output);
}

//------------------------------------------------
// reshuffle a vector
void shuffle1(vector<int> &x) {
    int rand1, dummy;
    int xSize = int(x.size());
    for (int i=0; i<xSize; i++) {
        rand1 = sample2(i,xSize-1);
        dummy = x[i];
        x[i] = x[rand1];
        x[rand1] = dummy;
    }
}

//------------------------------------------------
// draw from Bernoulli(p) distribution
int rbernoulli1(double p) {
    bernoulli_distribution dist_bernoulli(p);
    return(dist_bernoulli(generator));
}

//------------------------------------------------
// draw number of unique groups under Chinese restaurant process
int rCRPgroups(int n, double alpha) {
    if (n<=1)
        return(n);
    double p;
    int ngroups = 1;
    for (int i=1; i<n; i++) {
        p = alpha/double(i+alpha);
        if (rbernoulli1(p))
            ngroups++;
    }
    return(ngroups);
}

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd) {
    normal_distribution<double> rnorm(mean,sd);
    return(rnorm(generator));
}
