/* $COPYRIGHT_UM
   $LICENSE_LGPL
*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

#include "aprilsam.h"

//graph utils
void april_graph_factor_eval_destroy(april_graph_factor_eval_t *eval)
{
    if (!eval)
        return;

    for (int i = 0; i < eval->length; i++) {
        if (!eval->jacobians[i])
            break;

        matd_destroy(eval->jacobians[i]);
    }

    free(eval->jacobians);
    free(eval->r);
    matd_destroy(eval->W);
    free(eval);
}

int april_graph_dof(april_graph_t *graph)
{
    int factor_dof = 0, state_dof = 0;

    for (int i = 0; i < zarray_size(graph->nodes); i++) {
        april_graph_node_t *n;
        zarray_get(graph->nodes, i, &n);

        state_dof += n->length;
    }

    for (int i = 0; i < zarray_size(graph->factors); i++) {
        april_graph_factor_t *factor;
        zarray_get(graph->factors, i, &factor);

        factor_dof += factor->length;
    }

    return factor_dof - state_dof;
}

double april_graph_chi2(april_graph_t *graph)
{
    double chi2 = 0;

    for (int i = 0; i < zarray_size(graph->factors); i++) {
        april_graph_factor_t *factor;
        zarray_get(graph->factors, i, &factor);
        if(factor->type == APRIL_GRAPH_FACTOR_XYT_TYPE) {
            april_graph_factor_eval_t *eval = factor->state_eval(factor, graph, NULL);
            chi2 += 0.5 * eval->chi2;
            april_graph_factor_eval_destroy(eval);
        } else {
            april_graph_factor_eval_t *eval = factor->eval(factor, graph, NULL);
            chi2 += eval->chi2;
            april_graph_factor_eval_destroy(eval);
        }
    }

    return chi2;
}

// attributes
static void attr_add(april_graph_attr_t **attr, const stype_t *stype, const char *key, void *value)
{
    // XXX race
    if (*attr == NULL)
        *attr = april_graph_attr_create();

    zhash_t *hash = (*attr)->hash;

    struct april_graph_attr_record record = { .stype = stype,
                                              .value = value };
    char *keycopy = strdup(key);

    char *oldkey;
    struct april_graph_attr_record oldrecord;

    if (zhash_put(hash, &keycopy, &record, &oldkey, &oldrecord)) {
        // free the old value
        if (oldrecord.value == NULL) {
            // nothing to do when adding a value when only a NULL
            // exited before.
        } else if (oldrecord.value != value) {
            // they've changed the value. This application is
            // responsible for correctly removing the old one.
        }

        // free the old key
        free(oldkey);
    }
}

static void *attr_get(april_graph_attr_t *attr, const char *key)
{
    if (attr == NULL)
        return NULL;

    struct april_graph_attr_record record;
    if (zhash_get(attr->hash, &key, &record))
        return record.value;

    return NULL;
}

void april_graph_attr_put(april_graph_t *graph, const stype_t *type, const char *key, void *data)
{
    attr_add(&graph->attr, type, key, data);
}

void* april_graph_attr_get(april_graph_t * graph,
                           const char * key)
{
    return attr_get(graph->attr, key);
}

void april_graph_factor_attr_put(april_graph_factor_t *factor, const stype_t *type,
                                 const char *key, void *data)
{
    attr_add(&factor->attr, type, key, data);
}

void* april_graph_factor_attr_get(april_graph_factor_t * factor,
                                  const char * key)
{
    return attr_get(factor->attr, key);
}

void april_graph_node_attr_put(april_graph_node_t *node, const stype_t *type,
                               const char *key, void *data)
{
    attr_add(&node->attr, type, key, data);
}

void* april_graph_node_attr_get(april_graph_node_t *node,
                                const char *key)
{
    return attr_get(node->attr, key);
}

static void april_graph_attr_encode(const stype_t *stype, uint8_t *data, uint32_t *datapos, const void *obj)
{
    const april_graph_attr_t *attr = obj;

    zhash_iterator_t zit;
    zhash_iterator_init(attr->hash, &zit);
    char *name;
    struct april_graph_attr_record record;
    while (zhash_iterator_next(&zit, &name, &record)) {
        if (record.stype) {
            encode_u8(data, datapos, 1); // pompeil flag
            encode_string_u32(data, datapos, name);
            stype_encode_object(data, datapos, record.stype, record.value);
        } else {
//            printf("graph attribute %s doesn't have a type\n", name);
        }
    }

    encode_u8(data, datapos, 0);
}

static void *april_graph_attr_decode(const stype_t *stype, const uint8_t *data, uint32_t *datapos, uint32_t datalen)
{
    april_graph_attr_t *attr = april_graph_attr_create();

    while (decode_u8(data, datapos, datalen)) {
        struct april_graph_attr_record record;
        char *key = decode_string_u32(data, datapos, datalen);
        record.value = stype_decode_object(data, datapos, datalen, &record.stype);

        zhash_put(attr->hash, &key, &record, NULL, NULL);
    }

    return attr;
}

const stype_t stype_april_graph_attr = { .name = "april_graph_attr_t",
                                         .encode = april_graph_attr_encode,
                                         .decode = april_graph_attr_decode,
                                         .copy = NULL };

april_graph_attr_t *april_graph_attr_create()
{
    april_graph_attr_t *attr = calloc(1, sizeof(april_graph_attr_t));
    attr->stype = &stype_april_graph_attr;
    attr->hash = zhash_create(sizeof(char*), sizeof(struct april_graph_attr_record),
                              zhash_str_hash, zhash_str_equals);
    return attr;
}

void april_graph_attr_destroy(april_graph_attr_t *attr)
{
    if (!attr)
        return;

    zhash_destroy(attr->hash);
    free(attr);
}

void april_graph_encode(const stype_t *stype, uint8_t *data, uint32_t *datapos, const void *_obj)
{
    const april_graph_t *graph = _obj;

    for (int i = 0; i < zarray_size(graph->nodes); i++) {
        april_graph_node_t *node;
        zarray_get(graph->nodes, i, &node);

        if (node->stype) {
            encode_u8(data, datapos, 1); // encoding a node
            stype_encode_object(data, datapos, node->stype, node);
        } else {
            printf("node without stype\n");
        }
    }

    for (int i = 0; i < zarray_size(graph->factors); i++) {
        april_graph_factor_t *factor;
        zarray_get(graph->factors, i, &factor);

        if (factor->stype) {
            encode_u8(data, datapos, 2); // encoding a factor
            stype_encode_object(data, datapos, factor->stype, factor);
        } else {
            printf("factor without stype\n");
        }
    }

    encode_u8(data, datapos, 0); // EOF

    april_graph_attr_t *attr = graph->attr;
    stype_encode_object(data, datapos, attr ? attr->stype : NULL, attr);
}

void *april_graph_decode(const stype_t *stype, const uint8_t *data, uint32_t *datapos, uint32_t datalen)
{
    april_graph_t *graph = april_graph_create();

    while (1) {
        uint8_t op = decode_u8(data, datapos, datalen);
        if (op == 0)
            break;

        switch (op) {

            case 1: { // node
                april_graph_node_t *node = stype_decode_object(data, datapos, datalen, NULL);
                if (node != NULL) {
                    zarray_add(graph->nodes, &node);
                }
                break;
            }

            case 2: { // factor
                april_graph_factor_t *factor = stype_decode_object(data, datapos, datalen, NULL);
                if (factor != NULL) {
                    zarray_add(graph->factors, &factor);
                }
                break;
            }

            default: {
                printf("bad opcode %d\n", op);
                assert(0);
                break;
            }
        }
    }

    graph->attr = stype_decode_object(data, datapos, datalen, NULL);
    return graph;
}

const stype_t stype_april_graph = { .name = "april_graph_t",
                                    .encode = april_graph_encode,
                                    .decode = april_graph_decode,
                                    .copy = NULL };


april_graph_t *april_graph_create()
{
    april_graph_t *graph = calloc(1, sizeof(april_graph_t));

    graph->factors = zarray_create(sizeof(april_graph_factor_t*));
    graph->nodes = zarray_create(sizeof(april_graph_node_t*));
    graph->stype = &stype_april_graph;

    return graph;
}

void april_graph_destroy(april_graph_t *graph)
{
    for (int i = 0; i < zarray_size(graph->nodes); i++) {
        april_graph_node_t *n;
        zarray_get(graph->nodes, i, &n);

        n->destroy(n);
    }

    for (int i = 0; i < zarray_size(graph->factors); i++) {
        april_graph_factor_t *factor;
        zarray_get(graph->factors, i, &factor);

        factor->destroy(factor);
    }

    zarray_destroy(graph->nodes);
    zarray_destroy(graph->factors);

    free(graph);
}

void april_graph_xyt_stype_init();
void april_graph_xytpos_stype_init();

void april_graph_stype_init()
{
    stype_register(&stype_april_graph);
    stype_register(&stype_april_graph_attr);
    april_graph_xyt_stype_init();
    april_graph_xytpos_stype_init();
}

int april_graph_save(april_graph_t *graph, const char *path)
{
    uint32_t datalen = 0;
    stype_encode_object(NULL, &datalen, graph->stype, graph);

    uint8_t *tmp = malloc(datalen);
    uint32_t tmplen = 0;
    stype_encode_object(tmp, &tmplen, graph->stype, graph);

    FILE *f = fopen(path, "w");
    if(f == NULL)
    {
        printf("failed to open %s\n", path);
        return(0);
    }
    size_t res = fwrite(tmp, 1, tmplen, f);
    free(tmp);
    fclose(f);
    if (res == tmplen)
        return 1;
    return 0;
}

april_graph_t *april_graph_create_from_file(const char *path)
{
    april_graph_t *graph = NULL;

    FILE *f = fopen(path, "r");
    if (f == NULL)
        return NULL;

    fseek(f, 0L, SEEK_END);
    long len = ftell(f);
    fseek(f, 0L, SEEK_SET);

    uint8_t *tmp = malloc(len);

    size_t sz = fread(tmp, 1, len, f);
    if (sz != len)
        goto cleanup;

    uint32_t pos = 0;
    graph = stype_decode_object(tmp, &pos, len, NULL);

  cleanup:

    free(tmp);
    fclose(f);
    return graph;
}
