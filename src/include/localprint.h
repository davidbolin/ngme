//
//  localprint.h
//  mixture_library
//
//  Created by Jonas Wallin on 17/04/14.
//  Copyright (c) 2014 Jonas Wallin. All rights reserved.
//

#ifndef mixture_library_localprint_h
#define mixture_library_localprint_h

#ifdef matlab
	#include <mex.h>
	#define localprint(...) mexPrintf(__VA_ARGS__)
#else
	#define localprint(...) printf( __VA_ARGS__)
#endif

#endif
