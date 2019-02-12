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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "zarray.h"

int zstrcmp(const void * a_pp, const void * b_pp)
{
    assert(a_pp != NULL);
    assert(b_pp != NULL);

    char * a = *(void**)a_pp;
    char * b = *(void**)b_pp;

    return strcmp(a,b);
}

void zarray_vmap(zarray_t *za, void (*f)())
{
    assert(za != NULL);
    assert(f != NULL);
    assert(za->el_sz == sizeof(void*));

    for (int idx = 0; idx < za->size; idx++) {
        void *pp = &za->data[idx*za->el_sz];
        void *p = *(void**) pp;
        f(p);
    }
}
