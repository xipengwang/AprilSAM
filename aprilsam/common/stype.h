/* $COPYRIGHT_UM
$LICENSE_LGPL
*/

#ifndef _STYPE_H
#define _STYPE_H

#include <stdint.h>
#include <string.h>

#include "encode_bytes.h"

typedef struct stype stype_t;
struct stype
{
    char *name; // name of the type

    // encode and decode methods do not need to encode their own type
    // in any way.

    // Note that datapos can be NULL, which means that a "dry-run" is
    // performed. datapos is still updated, allowing a caller to know
    // how much data should be allocated.
    void (*encode)(const stype_t *stype, uint8_t *data, uint32_t *datapos, const void *obj);

    // During decode, an object does NOT read its type; it was already
    // consumed in order to find the right deserializer.
    void *(*decode)(const stype_t *stype, const uint8_t *data, uint32_t *datapos, uint32_t datalen);

    void *(*copy)(const stype_t *stype, const void *obj);

    void (*destroy)(const stype_t *stype, void *obj);

    void *impl;
};

// defines uint64 and friends
void stype_register_basic_types();

// the object "stype" should be allocated by caller and never modified
// or freed.
void stype_register(const stype_t *stype);

// returns NULL or the stype associated with name.
stype_t *stype_get(char *name);

void stype_encode_object(uint8_t *data, uint32_t *datapos, const stype_t *stype, const void *obj);

// outstype can be NULL.
void *stype_decode_object(const uint8_t *data, uint32_t *datapos, uint32_t datalen, const stype_t **outstype);

// write a single stype (and its children) to a file. Returns zero on success.
int stype_write_file(const stype_t *stype, void *obj, const char *path);

// read a single stype from a file, returns NULL on error.
void *stype_read_file(const char *path);

#endif
