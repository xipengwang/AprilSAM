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
$LICENSE_LGPL
*/

#include <stdio.h>
#include <stdint.h>
#include "stype.h"
#include "encode_bytes.h"

static void uint64_encode(const stype_t *stype, uint8_t *data, uint32_t *datapos, const void *obj)
{
    const uint64_t *p = obj;
    encode_u64(data, datapos, *p);
}

static void* uint64_decode(const stype_t *stype, const uint8_t *data, uint32_t *datapos, uint32_t datalen)
{
    uint64_t *p = malloc(sizeof(uint64_t));
    *p = decode_u64(data, datapos, datalen);
    return p;
}

static void *uint64_copy(const stype_t *stype, const void *obj)
{
    const int64_t *in = obj;
    int64_t *out = malloc(sizeof(uint64_t));
    *out = *in;
    return out;
}

static void uint64_destroy(const stype_t *stype, void *obj)
{
    free(obj);
}

const stype_t stype_uint64 = { .name = "uint64" ,
                              .encode = uint64_encode,
                              .decode = uint64_decode,
                              .copy = uint64_copy,
                              .destroy = uint64_destroy };

static void string_encode(const stype_t *stype, uint8_t *data, uint32_t *datapos, const void *obj)
{
    const char *p = obj;
    encode_string_u32(data, datapos, p);
}

static void *string_decode(const stype_t *stype, const uint8_t *data, uint32_t *datapos, uint32_t datalen)
{
    char *p;
    p = decode_string_u32(data, datapos, datalen);
    return p;
}

static void string_destroy(const stype_t *stype, void *obj)
{
    free(obj);
}

const stype_t stype_string = { .name = "string" ,
                               .encode = string_encode,
                               .decode = string_decode,
                               .copy = NULL,
                               .destroy = string_destroy };

void stype_register_basic_types()
{
    stype_register(&stype_uint64);
    stype_register(&stype_string);
}
