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

#ifndef _GETOPT_H
#define _GETOPT_H

#include "zarray.h"
#include "string_util.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct getopt getopt_t;

getopt_t *getopt_create();
void getopt_destroy(getopt_t *gopt);

// Parse args. Returns 1 on success
int getopt_parse(getopt_t *gopt, int argc, char *argv[], int showErrors);
void getopt_do_usage(getopt_t *gopt);

// Returns a string containing the usage. Must be freed by caller
char * getopt_get_usage(getopt_t *gopt);

void getopt_add_spacer(getopt_t *gopt, const char *s);
void getopt_add_bool(getopt_t *gopt, char sopt, const char *lname, int def, const char *help);
void getopt_add_int(getopt_t *gopt, char sopt, const char *lname, const char *def, const char *help);
void getopt_add_string(getopt_t *gopt, char sopt, const char *lname, const char *def, const char *help);
void getopt_add_double(getopt_t *gopt, char sopt, const char *lname, const char *def, const char *help);

const char *getopt_get_string(getopt_t *gopt, const char *lname);
int getopt_get_int(getopt_t *getopt, const char *lname);
int getopt_get_bool(getopt_t *getopt, const char *lname);
double getopt_get_double(getopt_t *getopt, const char *lname);
int getopt_was_specified(getopt_t *gopt, const char *lname);
const zarray_t *getopt_get_extra_args(getopt_t *gopt);

#ifdef __cplusplus
}
#endif


#endif
