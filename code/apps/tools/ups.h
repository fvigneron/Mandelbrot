//
//  ups.h
//  Multipliers
//
//  Created by Nicolae Mihalache on 10/01/2024.
//

#ifndef ups_h
#define ups_h

#include <stdio.h>

#include "ntypes.h"

#define UPS_STATUS_UNAVAILABLE   0
#define UPS_STATUS_ONLINE        1
#define UPS_STATUS_BATTERY       2
#define UPS_STATUS_ALARM         4
#define UPS_STATUS_OFF           8
#define UPS_STATUS_FAULT        16
#define UPS_STATUS_LOW_BATTERY  32

#define UPS_DEFAULT_NAME   "eaton3s"

/// A delay in ms for polling the state of the UPS.
#define UPS_DEFAULT_POLL_DELAY_MS    250
#define UPS_MIN_POLL_DELAY_MS         25
#define UPS_MAX_POLL_DELAY_MS       2500
#define UPS_NANO_SLEEP         250000000

bool ups_set_name(char *name);
bool ups_set_poll_delay_ms(ulong ms);

int ups_get_status(void);

void ups_wait_online(void);

#endif /* ups_h */
