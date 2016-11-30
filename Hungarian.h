//
//  MavericK
//  Hungarian.h
//
//  Created: Bob on 25/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Implements the Hungarian algorithm for finding the optimal assignment given a linear sum assignment problem. In this application this method is used to update label permutations on-the-fly to ensure a consistent labelling of components using the method of Stephens (2000).
//
//  Stephens, Matthew. "Dealing with label switching in mixture models." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 62.4 (2000) 795-809.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__Hungarian__
#define __Maverick1_0__Hungarian__

#include <iostream>
#include <random>
#include "readIn.h"

//------------------------------------------------
// the functions augmentLeft and augmentRight work together to find an augmented path. They call each other, which normally could lead to an infinite recursion, however, this is avoided as eventually either an augmented path will be found or no more moves will be possible. Return full path, or -1 if no path found.
std::vector<int> augmentLeft(int i, std::vector< std::vector<double> > &M, std::vector<int> &edgesRight, std::vector<int> &blockedLeft, std::vector<int> &blockedRight);
std::vector<int> augmentRight(int j, std::vector< std::vector<double> > &M, std::vector<int> &edgesRight, std::vector<int> &blockedLeft, std::vector<int> &blockedRight);

//------------------------------------------------
// carry out Hungarian algorithm to find best matching given cost matrix M
std::vector<int> hungarian(std::vector< std::vector<double> > &M, std::vector<int> &edgesLeft, std::vector<int> &edgesRight, std::vector<int> &blockedLeft, std::vector<int> &blockedRight);

#endif
