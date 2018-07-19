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
