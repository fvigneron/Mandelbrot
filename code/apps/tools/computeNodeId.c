//
//  computeNodeId.c
//  Mandelbrot
//
//  Created by Francois Vigneron on 23/12/2020.
//  Copyright © 2020 UPEC. All rights reserved.
//
#define _GNU_SOURCE  // must precede any other include if sched.h is included too
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "computeNodeId.h"

int get_cpu_id(void) {

#if defined(__linux__)
#include <sched.h>

    return sched_getcpu();
    
#else

    // Should run on linux kernel above 2.2.8 (if /proc/self/stat access is enabled)
    // For other arch, cpu_id defaults to -1
    long to_read = 8192;
    char buffer[to_read];
    int ok;
    int cpu_id = -1;
    
    // Get the the current process' stat file from the proc filesystem
    FILE* procfile = fopen("/proc/self/stat", "r");
    if( procfile != NULL ) {
        ok = fread(buffer, sizeof(char), to_read, procfile) == to_read;
    } else {
        ok = 0;
    }
    fclose(procfile);
    
    // Looking for the 39th entry
    if ( ok ) { char* line = strtok(buffer, " ");
        int i;
        for (i = 1; i < 38; i++) {
            line = strtok(NULL, " ");
        }
        line = strtok(NULL, " ");
        cpu_id = atoi(line);
    }
    
    return cpu_id;
    
#endif

}
