//
//  MavericK
//  Windows_functions.h
//
//  Created: Bob on 18/05/2016
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Contains functions that are present by default when compiling on Xcode 7.3.1, but absent from Visual Studio 2010.
//  Code for gamma() and lgamma() functions taken from original code due to John D. Cook, available at http://www.johndcook.com/stand_alone_code.html and distributed under BSD license.
//  Many thanks John!
//  Note that the functions gamma() and lgamma() are mutually dependent.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__Windows_functions__
#define __Maverick1_0__Windows_functions__

#ifdef _WIN32
//------------------------------------------------
// round a double to nearest integer
int round(double x);

//------------------------------------------------
// the gamma function
double gamma(double);

//------------------------------------------------
// the natural logarithm of the gamma function
double lgamma(double);

#endif
#endif

