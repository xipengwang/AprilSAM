/* $COPYRIGHT_UM
$LICENSE_LGPL
*/

#ifndef _ZSET_H
#define _ZSET_H

#include "zhash.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef zhash_t zset_t;
typedef zhash_iterator_t zset_iterator_t;

static inline zset_t *zset_create(size_t keysz,
                                   uint32_t(*hash)(const void *a),
                                   int(*equals)(const void *a, const void *b))
{
    return zhash_create(keysz, 0, hash, equals);
}

static inline zset_t *zset_copy(zset_t *zs)
{
    return zhash_copy(zs);
}

static inline void zset_vmap(zset_t *zs, void (*f)())
{
    zhash_vmap_keys(zs, f);
}

static inline void zset_map(zset_t *zs, void (*f)())
{
    zhash_map_keys(zs, f);
}

static inline void zset_destroy(zset_t *zs)
{
    zhash_destroy(zs);
}

static inline int zset_add(zset_t *zs, const void *key, void *oldkey)
{
    const void *value = NULL; // with valuesz zero, value won't be touched
    return zhash_put(zs, key, &value, oldkey, NULL);
}

static inline int zset_remove(zset_t *zs, const void *key, void *oldkey)
{
    return zhash_remove(zs, key, oldkey, NULL);
}

static inline void zset_clear(zset_t *zs)
{
    zhash_clear(zs);
}

static inline int zset_contains(zset_t *zs, const void *key)
{
    return zhash_contains(zs, key);
}

static inline int zset_size(zset_t *zs)
{
    return zhash_size(zs);
}

// do not modify the set while iterating over the elements.
static inline void zset_iterator_init(zset_t *zs, zset_iterator_t *zit)
{
    zhash_iterator_init(zs, zit);
}

static inline int zset_iterator_next(zset_iterator_t *zit, void *outkey)
{
    return zhash_iterator_next(zit, outkey, NULL);
}

/** Returns non-zero if each set is the same size and contains the
 *  same elements. 'scratch_ptr' points to memory of size 'keysz' so
 *  zset_equals has a place to temporarily store keys
 */
static inline int zset_equals(zset_t *a, zset_t *b, void *scratch_ptr)
{
    zset_iterator_t it;

    if (zset_size(a) != zset_size(b))
        return 0;

    zset_iterator_init(a, &it);
    while (zset_iterator_next(&it, scratch_ptr)) {
        if (!zset_contains(b, scratch_ptr))
            return 0;
    }

    zset_iterator_init(b, &it);
    while (zset_iterator_next(&it, scratch_ptr)) {
        if (!zset_contains(a, scratch_ptr))
            return 0;
    }

    return 1;
}

#ifdef __cplusplus
}
#endif

#endif
