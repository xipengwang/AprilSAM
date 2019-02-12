/* Copyright (C) 2013-2019, The Regents of The University of Michigan.
All rights reserved.

This software was developed in the APRIL Robotics Lab under the
direction of Edwin Olson, ebolson@umich.edu. This software may be
available under alternative licensing terms; contact the address above.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, see <http://www.gnu.org/licenses/>.
*/

/*$LICENSE*/

/* $COPYRIGHT_UM
$LICENSE_BSD
*/

#ifndef _TIME_UTIL_H
#define _TIME_UTIL_H

#include <stdbool.h>
#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct timeutil_rest timeutil_rest_t;
timeutil_rest_t *timeutil_rest_create();
void timeutil_rest_destroy(timeutil_rest_t * rest);

int64_t utime_now(); // blacklist-ignore
int64_t utime_get_seconds(int64_t v);
int64_t utime_get_useconds(int64_t v);
void    utime_to_timeval(int64_t v, struct timeval *tv);
void    utime_to_timespec(int64_t v, struct timespec *ts);

int32_t  timeutil_usleep(int64_t useconds);
uint32_t timeutil_sleep(unsigned int seconds);
int32_t  timeutil_sleep_hz(timeutil_rest_t *rest, double hz);

void timeutil_timer_reset(timeutil_rest_t *rest);
void timeutil_timer_start(timeutil_rest_t *rest);
void timeutil_timer_stop(timeutil_rest_t *rest);
bool timeutil_timer_timeout(timeutil_rest_t *rest, double timeout_s);

int64_t time_util_hhmmss_ss_to_utime(double time);



#ifdef __cplusplus
}
#endif

#endif
