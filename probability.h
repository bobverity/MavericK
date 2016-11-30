//
//  MavericK
//  probability.h
//
//  Created: Bob on 22/09/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Functions for sampling from some fairly basic probability mass and density functions.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__probability__
#define __Maverick1_0__probability__

#include <random>
#include "misc.h"

//------------------------------------------------
// draw from uniform(a,b) distribution
double runif1(double a=0, double b=1);

//------------------------------------------------
// draw from gamma(shape,rate) distribution
double rgamma1(double shape, double rate);

//------------------------------------------------
// draw from beta(shape1,shape2) distribution
double rbeta1(double shape1, double shape2);

//------------------------------------------------
// sample from given probability vector (that sums to pSum)
int sample1(std::vector<double> &p, double pSum);

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal probability
int sample2(int b, int a=1);

//------------------------------------------------
// sample n terms without replacement from vector x
std::vector<int> sample3(std::vector<int> x, int n);

//------------------------------------------------
// reshuffle a vector
void shuffle1(std::vector<int> &x);

//------------------------------------------------
// draw from Bernoulli(p) distribution
int rbernoulli1(double p);

//------------------------------------------------
// draw number of unique groups under Chinese restaurant process
int rCRPgroups(int n, double alpha);

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd);

#endif
