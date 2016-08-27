//
//  TimeIt.h
//  mixture_library
//
//  Created by Jonas Wallin on 18/04/14.
//  Copyright (c) 2014 Jonas Wallin. All rights reserved.
//

#ifndef mixture_library_TimeIt_h
#define mixture_library_TimeIt_h

#include <ctime>


#define timer_start(num)  const double ticks_per_ms_##num = static_cast<double>(CLOCKS_PER_SEC); \
                          clock_t start_##num = clock(); 
#define timer_end(num,...)    double time_init_##num = static_cast<double>(clock()-start_##num)  / ticks_per_ms_##num; \
                            printf("%s =  %.8f seconds.\n",__VA_ARGS__, time_init_##num);

#define timer_init_loop(num) double START_TIMERS__[num] = { 0 }; const double TICKS_PER_MS__ = static_cast<double>(CLOCKS_PER_SEC);
#define timer_loop_start(num)  START_TIMERS__[num-1] -= clock();
#define timer_loop_end(num)  START_TIMERS__[num-1] += clock();
#define timer_print(num,...)  printf("%s =  %.8f seconds.\n",__VA_ARGS__, START_TIMERS__[num-1]/ TICKS_PER_MS__)  ;


#endif
