//
//  MavericK
//  OSfunctions.h
//
//  Created: Bob on 22/09/2015
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Defines certain functions and commands that are OS-specific. Also conditionally reads in header files containing functions that are needed when compiling on particular operating systems.
//
// ---------------------------------------------------------------------------

#ifndef __Maverick1_0__OSfunctions__
#define __Maverick1_0__OSfunctions__

#ifdef _WIN32 // if Windows
    #include <direct.h>
	#define my_isnan _isnan
	#define my_gamma approx_gamma
	#define my_lgamma approx_lgamma
    #define DIRBREAK "\\"
    #define GETCWD _getcwd
	#define isnan(x) _isnan(x)

#else // if not Windows
    #include <unistd.h>
	#define my_isnan std::isnan
	#define my_gamma std::gamma
	#define my_lgamma std::lgamma
    #define DIRBREAK "/"
    #define GETCWD getcwd
#endif

#endif
