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
// sample from given probability vector (that sums to pSum)
int sample1(std::vector<double> &p, double pSum);

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd);

#endif
