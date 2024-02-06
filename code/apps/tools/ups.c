//
//  ups.c
//  Multipliers
//
//  Created by Nicolae Mihalache on 10/01/2024.
//

#include <string.h>

#include "ntypes.h"
#include "ups.h"
#include "stopWatch.h"
#include "io.h"

static volatile bool ups_available  = true, def_cmd = false;
static char ups_cmd[200];
static volatile ulong ups_delay = UPS_DEFAULT_POLL_DELAY_MS;
static volatile ulong ups_ts = 0;

static bool all_caps(char *line) {
    for (int i = 0; line[i] != 0; i++) {
        char c = line[i];
        if(c == ' ' || c == '\n') {
            continue;
        }
        
        if(c < 'A' || c > 'Z') {
            return false;
        }
    }
    
    return true;
}

bool ups_set_name(char *name) {
    if(name == NULL || strlen(name) > 150) {
        return false;
    }
    
    snprintf(ups_cmd, 200, "upsc %s ups.status", name);
    ups_available = true;
    def_cmd = true;
    
    return true;
}

bool ups_set_poll_delay_ms(ulong ms) {
    if(ms < UPS_MIN_POLL_DELAY_MS || ms > UPS_MAX_POLL_DELAY_MS) {
        return false;
    }
    
    ups_delay = ms;
    
    return true;
}

int ups_get_status(void) {
    if(! def_cmd) {
        ups_set_name(UPS_DEFAULT_NAME);
    }
    
    if(! ups_available) {
        return UPS_STATUS_UNAVAILABLE;
    }
    
    FILE *pipe;
    char line[200];
    
    pipe = popen(ups_cmd, "r");
    if(pipe == NULL || fgets(line, 200, pipe) == NULL || ! all_caps(line)) {
        ups_available = false;
        
        return UPS_STATUS_UNAVAILABLE;
    }
    pclose(pipe);
    
    int status = 0;
    
    if(strstr(line, "OL") != NULL) {
        status |= UPS_STATUS_ONLINE;
    }
    
    if(strstr(line, "OB") != NULL) {
        status |= UPS_STATUS_BATTERY;
    }
    
    if(strstr(line, "ALARM") != NULL) {
        status |= UPS_STATUS_ALARM;
    }
    
    if(strstr(line, "OFF") != NULL) {
        status |= UPS_STATUS_OFF;
    }
    
    if(strstr(line, "RB") != NULL) {
        status |= UPS_STATUS_FAULT;
    }
    
    if(strstr(line, "LB") != NULL) {
        status |= UPS_STATUS_LOW_BATTERY;
    }
    
    ups_available = status != 0;
    
    return status;
}

void ups_wait_online(void) {
    ptime now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    ulong ts = now.tv_sec * 1000 + now.tv_nsec / 1000000;
    if(ts - ups_ts >= ups_delay) {
        ups_ts = ts;
        
        char cdat[100];
        int ups = ups_get_status();
        bool delay = false;
        
        if(ups != UPS_STATUS_UNAVAILABLE && ups != UPS_STATUS_ONLINE) {
            date(cdat, 100);
            io_log_echo("UPS status: %d, delaying threads (%s)\n", ups, cdat);
            delay = true;
        }
        
        struct timespec sleep = {.tv_sec = 0, .tv_nsec = UPS_NANO_SLEEP};
        while(ups != UPS_STATUS_UNAVAILABLE && ups != UPS_STATUS_ONLINE) {
            nanosleep(&sleep, NULL);
            ups = ups_get_status();
        }
        
        if(delay) {
            date(cdat, 100);
            io_log_echo("UPS back online, restarting threads (%s)\n", cdat);
        }
    }
}
