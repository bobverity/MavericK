//
//  MavericK
//  Notes.c
//
//  Created: Bob on 22/03/2016
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
/*
------------------------------------------------
OVERALL PROGRAM STRUCTURE
 
- The overall flow of the program is specified in main.cpp. First, global parameters are defined (given default values) in a single "globals" object. Then a handful of input arguments are processed - just those needed to locate the parameters.txt file. This file is then located and read in, and finally remaining command-line arguments are processed. Therefore the order of precedence is: default values, which are overwritten by parameters in the parameters.txt file, which are overwritten by parameters specified on the command line.
- The data.txt file is then read in, and certain checks are carried out to ensure that the combination of data and parameters makes sense.
- The main loop of the program over values of K then begins. In each iteration there are various different sorts of analysis carried out (for example exhaustive analysis, thermodynamic integration etc.). The main MCMC loop is compulsory, but all other types of analysis are optional.
- For MCMC analyses (either the ordinary MCMC or the thermodynamic integral MCMC, and either with- or without-admixture) a distinct MCMCobject class is defined. This class contains all functions and variables needed to carry out the MCMC. Final results are saved back into the globals object.
- Where possible, outputs are produced at each iteration of the main loop (for each K) and results are flushed. This means that outputs can be viewed even when the program has not completely finished. Some outputs (such as normalised evidence values) can only be produced after the main loop.
- This code is designed to compile on both Windows and Mac. Operating-system specific functions are contained in the header file "OSfunctions.h"

------------------------------------------------
EXHAUSTIVE APPROACH - SUMMING OVER PARTITIONS ACCORDING TO RESTRICTED GROWTH FUNCTION

 The exact model evidence can be obtained by summing over all possible assignments of n objects to K groups (n is the number of individuals in the case of the without-admixture model, or the number of gene copies in the case of the admixture model). For example, if n=6 and K=3 then we would need to explore the following allocations...
 
 1 1 1 1 1 1
 1 1 1 1 1 2
 1 1 1 1 1 3
 1 1 1 1 2 1
 1 1 1 1 2 2
 1 1 1 1 2 3
 1 1 1 1 3 1
 ...
 3 3 3 3 3 2
 3 3 3 3 3 3
 
 The number of such allocations is K^n, which becomes prohibitively large for even moderate n and K.
 We can make this more efficient by only summing over unique partitions. Mixture components are not uniquely defined, meaning allocations that correspond to the same partition have the same likelihood. For example, {1 1 1 2 2 2} is equivalent to {2 2 2 1 1 1} is equivalent to {3 3 3 2 2 2} etc. We can loop over all unique partitions using the restricted growth function, in which the ith element can never be more than 1 higher than preceeding values (and no value can exceed K). This leads to the following list of allocations...
 
 1 1 1 1 1 1
 1 1 1 1 1 2
 1 1 1 1 2 1
 1 1 1 1 2 2
 1 1 1 1 2 3
 1 1 1 2 1 1
 1 1 1 2 1 2
 1 1 1 2 1 3
 1 1 1 2 2 1
 1 1 1 2 2 2
 1 1 1 2 2 3
 ...
 1 2 3 3 3 2
 1 2 3 3 3 3
 
 The number of such partitions is given by the sum from k=1 to k=K of S{n,k}, where S{n,k} denotes the Stirling numbers of the second kind (note that when K>=n this is equivalent to the Bell number B_n). This is often significantly smaller than K^n.
When calculating the overall likelihood, the likelihood of each partition must be multiplied by the number of unique group allocations that it encapsulates. For example, the partition {1 1 1 1 1 1} encapsulates K unique allocations, because the "1" values could be replaced with any value from 1 to K and still correspond to this partition.
Taking a more complex example, the partition {1 1 1 2 2 3} can be written as follows: {o o o}{o o}{o}, where each dot "o" can be any unique number. There are K values that we can assign to the first cluster of three dots, then (K-1) values that we can assign to the second cluster of two dots, then (K-2) values that we can assign to the final cluster of a single dot. So, the number of unique allocations represented by this partition is K(K-1)(K-2).
 More generally, let u denote the total number of elements in the partition (also the largest value in the restricted growth function). So in the above example u=3. Then the number of allocation is given by K!/(K-u)!.
This is why log-likelihoods are weighted by lgamma(double(K)+1)-lgamma(double(K)-double(uniques)+1) in the exhaustive approach.


 ------------------------------------------------
 MIT License
 
 Copyright (c) 2016 Robert Verity
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*/