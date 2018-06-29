/* $COPYRIGHT_UM
$LICENSE_LGPL
*/

#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include "stype.h"
#include "zset.h"
#include "io_util.h"

#define TNAME stype_hash
#define TKEYTYPE char*
#define TVALTYPE stype_t*

// takes a pointer to the key
#define TKEYHASH(pk) ((*pk)[0])
// takes a pointer to the value
#define TKEYEQUAL(pka, pkb) (!strcmp(*pka, *pkb))
#include "thash_impl.h"
#undef TNAME
#undef TKEYTYPE
#undef TVALTYPE

static stype_hash_t *stypes;

void stype_register(const stype_t *stype)
{
    // XXX race
    if (stypes == NULL) {
        stypes = stype_hash_create();
    }

    stype_hash_put(stypes, (char**) &stype->name, (stype_t**) &stype, NULL, NULL);
}

stype_t *stype_get(char *name)
{
    stype_t *stype = NULL;

    if (stypes == NULL)
        return NULL;

    if (stype_hash_get(stypes, &name, &stype))
        return stype;

    return NULL;
}

// you CAN encode NULL (in which case stype MAY be null also)
void stype_encode_object(uint8_t *data, uint32_t *datapos, const stype_t *stype, const void *obj)
{
    assert(obj == NULL || stype != NULL);

    // XXX eventually, randomly generate these cookies?
    static uint64_t next_magic = 0x7b287f8a1579a0edULL;

    uint64_t magic = next_magic;
    next_magic++;

    encode_u64(data, datapos, magic);

    if (obj == NULL) {
        encode_string_u32(data, datapos, "");
        encode_u32(data, datapos, 0);
    } else {
        assert(stype->name != NULL);

        encode_string_u32(data, datapos, stype->name);

        // compute length and serialize
        uint32_t datapos0 = *datapos;
        stype->encode(stype, NULL, datapos, obj);
        uint32_t length = *datapos - datapos0;
        *datapos = datapos0;

        encode_u32(data, datapos, length);

        stype->encode(stype, data, datapos, obj);
    }

    encode_u64(data, datapos, magic);
}

void *stype_decode_object(const uint8_t *data, uint32_t *datapos, uint32_t datalen, const stype_t **outstype)
{
    static zset_t *no_stype_warnings = NULL;

    uint64_t magic = decode_u64(data, datapos, datalen);
    char *stype_name = decode_string_u32(data, datapos, datalen);
    uint32_t length = decode_u32(data, datapos, datalen);

    void *obj = NULL;

    if (1) { // used to be length > 0, but a zero-length object is okay.
        stype_t *stype = stype_get(stype_name);

        if (stype != NULL) {
            obj = stype->decode(stype, data, datapos, datalen);
            if (outstype)
                *outstype = stype;

        } else if (length > 0) {
            // TODO seek until we find another copy of 'magic'
            if (!no_stype_warnings) {
                no_stype_warnings = zset_create(sizeof(char*),
                                                zhash_str_hash, zhash_str_equals);
            }
            if (!zset_contains(no_stype_warnings, &stype_name)) {
                char *s = strdup(stype_name);
                zset_add(no_stype_warnings, &s, NULL);
                printf("Unknown stype %s\n", stype_name);
            }

            if (outstype)
                *outstype = NULL;
            int okay = 0;

            for (uint32_t pos = *datapos; pos + 8 <= datalen; pos++) {
                uint32_t tmp = pos;
                uint32_t maybe_magic = decode_u64(data, &tmp, datalen);
                if (maybe_magic == magic) {
                    *datapos = pos;
                    okay = 1;
                    break;
                }
            }

            if (!okay) {
                printf("failed to seek to next object.\n");
            }
        }
    }

    uint64_t magic2 = decode_u64(data, datapos, datalen);

    if (magic != magic2) {
        printf("magic mismatch while decoding '%s': %016"PRIx64" %016"PRIx64"\n", stype_name, magic, magic2);
        assert(0);
    }

    free(stype_name);

    return obj;
}

int stype_write_file(const stype_t *stype, void *obj, const char *path)
{
    uint32_t len = 0;
    stype_encode_object(NULL, &len, stype, obj);

    uint8_t *buf = calloc(1, len);
    len = 0;
    stype_encode_object(buf, &len, stype, obj);

    int ret = ioutils_write_file(path, buf, len);
    free(buf);
    return ret;
}

void *stype_read_file(const char *path)
{
    char *buf;
    long bufsz;

    if (ioutils_read_file(path,  (void**) &buf, &bufsz, -1)) {
        return NULL;
    }

    uint32_t datapos = 0;
    void *obj = stype_decode_object((uint8_t*) buf, &datapos, bufsz, NULL);
    free(buf);
    return obj;
}
