//
//  MavericK
//  Hungarian.cpp
//
//  Created: Bob on 25/10/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include "Hungarian.h"

using namespace std;

//------------------------------------------------
// the functions augmentLeft and augmentRight work together to find an augmented path. They call each other, which normally could lead to an infinite recursion, but this is avoided as eventually either an augmented path will be found or no more moves will be possible. Return full path, or -1 if no path found.
vector<int> augmentLeft(int i, vector< vector<double> > &M, vector<int> &edgesRight, vector<int> &blockedLeft, vector<int> &blockedRight) {
    
    blockedLeft[i] = 1;
    vector<int> output(1,-1);
    
    // search all unmatched edges
    for (int j=0; j<int(M.size()); j++) {
        if (M[i][j]==0 && edgesRight[j]!=i && blockedRight[j]==0) {
            
            // if edge leads to augmented path then add current node to path and return
            output = augmentRight(j, M, edgesRight, blockedLeft, blockedRight);
            if (output[0]>=0) {
                output.push_back(i);
                return(output);
            }
        }
    }
    
    // if no more moves then return -1
    return(output);
}

vector<int> augmentRight(int j, vector< vector<double> > &M, vector<int> &edgesRight, vector<int> &blockedLeft, vector<int> &blockedRight) {
    
    blockedRight[j] = 1;
    vector<int> output(1);
    
    // if node j is unmatched then return j as start of augmented path
    if (edgesRight[j]<0) {
        output[0] = j;
        return(output);
    }
    
    // otherwise continue chain of augmenting
    output = augmentLeft(edgesRight[j], M, edgesRight, blockedLeft, blockedRight);
    if (output[0]>=0) {
        output.push_back(j);
    }
    return(output);
}

//------------------------------------------------
// carry out Hungarian algorithm to find best matching given cost matrix M
vector<int> hungarian(vector< vector<double> > &M, vector<int> &edgesLeft, vector<int> &edgesRight, vector<int> &blockedLeft, vector<int> &blockedRight) {
    int n = int(M.size());
    
    // define maximum number of reps in Hungarian algorithm before aborting
    int maxReps = int(1e6);
    
    // initialise assignment objects
    int numberAssigned;
    
    // search for solution until maxReps reached
    for (int rep=0; rep<maxReps; rep++) {
        
        // zero assignment objects
        for (int i=0; i<n; i++) {
            edgesLeft[i] = -1;
            edgesRight[i] = -1;
        }
        numberAssigned = 0;
        
        // subtract smallest element from all rows and columns
        double minRow;
		vector<double> minCol = M[0];
        for (int i=0; i<int(M.size()); i++) {
            minRow = *min_element(begin(M[i]),end(M[i]));
            for (int j=0; j<int(M[i].size()); j++) {
                M[i][j] -= minRow;
                if (M[i][j]<minCol[j])
                    minCol[j] = M[i][j];
            }
        }
        for (int i=0; i<int(M.size()); i++) {
            for (int j=0; j<int(M[i].size()); j++) {
                M[i][j] -= minCol[j];
            }
        }
        
        // generate an initial matching
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (M[i][j]==0 && edgesRight[j]<0) {
                    edgesLeft[i] = j;
                    edgesRight[j] = i;
                    numberAssigned++;
                    break;
                }
            }
        }
        
        // if this matching is perfect then we are done
        if (numberAssigned==n)
            return(edgesLeft);
        
        // continue augmenting paths until no more possible
        //vector<int> blockedLeft(n);
        //vector<int> blockedRight(n);
        bool continueAugmenting = true;
        while (continueAugmenting) {
            continueAugmenting = false;
            
            // search all unmatched nodes
            for (int i=0; i<n; i++) {
                if (edgesLeft[i]<0) {
                    
                    // attempt to find augmented path
                    blockedLeft = vector<int>(n);
                    blockedRight = vector<int>(n);
                    vector<int> path = augmentLeft(i, M, edgesRight, blockedLeft, blockedRight);
                    
                    // if successful then augment
                    if (path[0]>=0) {
                        continueAugmenting = true;
                        numberAssigned ++;
                        for (int j=0; j<int(path.size()/2); j++) {
                            edgesLeft[path[j*2+1]] = path[j*2];
                            edgesRight[path[j*2]] = path[j*2+1];
                        }
                        
                        // if best matching found then finish
                        if (numberAssigned==n) {
                            return(edgesLeft);
                        }
                    }
                }
            }
        }
        
        // find minimum value in cost matrix, looking at all elements in which neither the row or the column is part of the minimum vertex cover
        double minVal = -log(double(0));
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (blockedLeft[i]==1 && blockedRight[j]==0 && M[i][j]<minVal) {
                    minVal = M[i][j];
                }
            }
        }
        
        // add or subtract this value from cost matrix as required
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (blockedLeft[i]==1 && blockedRight[j]==0) {
                    M[i][j] -= minVal;
                }
                if (blockedLeft[i]==0 && blockedRight[j]==1) {
                    M[i][j] += minVal;
                }
            }
        }
        
        // at this point we have a new cost matrix and can repeat the process from the top
        
    } // rep loop
    
    exit(1);
}
