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

#include "aprilsam.h"
#include "common/doubles.h"
#include "common/matd.h"

/////////////////////////////////////////////////////////////////////////////////////////
// XYT Factor
static april_graph_factor_t* xyt_factor_copy(april_graph_factor_t *factor)
{
    april_graph_factor_t * next = calloc(1, sizeof(april_graph_factor_t));
    next->stype = factor->stype;
    next->type = factor->type;
    next->nnodes = factor->nnodes;
    next->nodes = calloc(1, sizeof(int) * next->nnodes);
    memcpy(next->nodes, factor->nodes, sizeof(int) * next->nnodes);
    next->length = factor->length;

    if (factor->attr) {
        //assert(0);
        next->attr = calloc(1, sizeof(april_graph_attr_t));
        next->attr->stype = factor->attr->stype;
        next->attr->hash = zhash_copy(factor->attr->hash);
    }
    next->copy = factor->copy;
    next->eval = factor->eval;
    next->state_eval = factor->state_eval;
    next->destroy = factor->destroy;

    next->u.common.z = doubles_dup(factor->u.common.z, 3);
    next->u.common.ztruth = doubles_dup(factor->u.common.ztruth, 3);
    next->u.common.W = matd_copy(factor->u.common.W);

    return next;
}

static april_graph_factor_eval_t* xyt_factor_eval(april_graph_factor_t *factor, april_graph_t *graph, april_graph_factor_eval_t *eval)
{
    if (eval == NULL) {
        eval = calloc(1, sizeof(april_graph_factor_eval_t));
        eval->jacobians = calloc(3, sizeof(matd_t*)); // NB: NULL-terminated
        eval->jacobians[0] = matd_create(3,3);
        eval->jacobians[1] = matd_create(3,3);
        eval->r = calloc(3, sizeof(double));
        eval->W = matd_create(3,3);
    }

    april_graph_node_t *na, *nb;
    zarray_get(graph->nodes, factor->nodes[0], &na);
    zarray_get(graph->nodes, factor->nodes[1], &nb);

    double xa = na->l_point[0], ya = na->l_point[1], ta = na->l_point[2];
    double xb = nb->l_point[0], yb = nb->l_point[1], tb = nb->l_point[2];

    double ca = cos(ta), sa = sin(ta);

    // predicted obs
    double dx = xb - xa, dy = yb - ya;
    double zhat[3];
    zhat[0] =  ca*dx + sa*dy;
    zhat[1] = -sa*dx + ca*dy;
    zhat[2] =  tb    - ta;

    // partial derivatives of zhat WRT node 0 [xa ya ta]
    matd_set_data(eval->jacobians[0], (double[]) {
               -ca,  -sa, -sa*dx+ca*dy,
                sa,  -ca, -ca*dx-sa*dy,
                0,     0,   -1 });

    // partial derivatives of zhat WRT node 1 [xb yb tb]
    // TODO: Shouldn't this be negative of J[0]? Why right column mostly zero?
    matd_set_data(eval->jacobians[1], (double[]) {
                ca,  sa, 0,
                -sa, ca, 0,
                0,    0, 1 });


    eval->length = 3;

    eval->r[0] = factor->u.common.z[0] - zhat[0];
    eval->r[1] = factor->u.common.z[1] - zhat[1];
    eval->r[2] = mod2pi(factor->u.common.z[2] - zhat[2]);

    memcpy(eval->W->data, factor->u.common.W->data, 3*3*sizeof(double));

    // chi^2 = r'*W*r, via X = W*r
    double X[3] = { MATD_EL(eval->W, 0, 0)*eval->r[0] +
                    MATD_EL(eval->W, 0, 1)*eval->r[1] +
                    MATD_EL(eval->W, 0, 2)*eval->r[2],
                    MATD_EL(eval->W, 1, 0)*eval->r[0] +
                    MATD_EL(eval->W, 1, 1)*eval->r[1] +
                    MATD_EL(eval->W, 1, 2)*eval->r[2],
                    MATD_EL(eval->W, 2, 0)*eval->r[0] +
                    MATD_EL(eval->W, 2, 1)*eval->r[1] +
                    MATD_EL(eval->W, 2, 2)*eval->r[2] };
    eval->chi2 = eval->r[0]*X[0] + eval->r[1]*X[1] + eval->r[2]*X[2];

    return eval;
}

static april_graph_factor_eval_t* xyt_factor_state_eval(april_graph_factor_t *factor, april_graph_t *graph, april_graph_factor_eval_t *eval)
{
    if (eval == NULL) {
        eval = calloc(1, sizeof(april_graph_factor_eval_t));
        eval->jacobians = calloc(3, sizeof(matd_t*)); // NB: NULL-terminated
        eval->jacobians[0] = matd_create(3,3);
        eval->jacobians[1] = matd_create(3,3);
        eval->r = calloc(3, sizeof(double));
        eval->W = matd_create(3,3);
    }

    april_graph_node_t *na, *nb;
    zarray_get(graph->nodes, factor->nodes[0], &na);
    zarray_get(graph->nodes, factor->nodes[1], &nb);

    double xa = na->state[0], ya = na->state[1], ta = na->state[2];
    double xb = nb->state[0], yb = nb->state[1], tb = nb->state[2];

    double ca = cos(ta), sa = sin(ta);

    // predicted obs
    double dx = xb - xa, dy = yb - ya;
    double zhat[3];
    zhat[0] =  ca*dx + sa*dy;
    zhat[1] = -sa*dx + ca*dy;
    zhat[2] =  tb    - ta;

    // partial derivatives of zhat WRT node 0 [xa ya ta]
    matd_set_data(eval->jacobians[0], (double[]) {
               -ca,  -sa, -sa*dx+ca*dy,
                sa,  -ca, -ca*dx-sa*dy,
                0,     0,   -1 });

    // partial derivatives of zhat WRT node 1 [xb yb tb]
    // TODO: Shouldn't this be negative of J[0]? Why right column mostly zero?
    matd_set_data(eval->jacobians[1], (double[]) {
                ca,  sa, 0,
                -sa, ca, 0,
                0,    0, 1 });


    eval->length = 3;

    eval->r[0] = factor->u.common.z[0] - zhat[0];
    eval->r[1] = factor->u.common.z[1] - zhat[1];
    eval->r[2] = mod2pi(factor->u.common.z[2] - zhat[2]);

    memcpy(eval->W->data, factor->u.common.W->data, 3*3*sizeof(double));

    // chi^2 = r'*W*r, via X = W*r
    double X[3] = { MATD_EL(eval->W, 0, 0)*eval->r[0] +
                    MATD_EL(eval->W, 0, 1)*eval->r[1] +
                    MATD_EL(eval->W, 0, 2)*eval->r[2],
                    MATD_EL(eval->W, 1, 0)*eval->r[0] +
                    MATD_EL(eval->W, 1, 1)*eval->r[1] +
                    MATD_EL(eval->W, 1, 2)*eval->r[2],
                    MATD_EL(eval->W, 2, 0)*eval->r[0] +
                    MATD_EL(eval->W, 2, 1)*eval->r[1] +
                    MATD_EL(eval->W, 2, 2)*eval->r[2] };
    eval->chi2 = eval->r[0]*X[0] + eval->r[1]*X[1] + eval->r[2]*X[2];

    return eval;
}

static void xyt_factor_destroy(april_graph_factor_t *factor)
{
    free(factor->nodes);
    free(factor->u.common.z);
    free(factor->u.common.ztruth);
    matd_destroy(factor->u.common.W);
    if(factor->attr) {
        zhash_iterator_t zit;
        april_graph_attr_t *attr = factor->attr;
        zhash_iterator_init(attr->hash, &zit);
        char *name;
        struct april_graph_attr_record record;
        while (zhash_iterator_next(&zit, &name, &record)) {
            if (record.stype) {
                record.stype->destroy(record.stype, record.value);
            } else {
                printf("graph attribute %s doesn't have a type\n", name);
            }
        }
        zhash_map_keys(attr->hash, free_key);
        zhash_destroy(attr->hash);
        free(factor->attr);
    }
    free(factor);
}

static void april_graph_factor_xyt_encode(const stype_t *stype, uint8_t *data, uint32_t *datapos, const void *obj)
{
    const april_graph_factor_t *factor = obj;

    encode_u32(data, datapos, factor->nodes[0]);
    encode_u32(data, datapos, factor->nodes[1]);

    for (int i = 0; i < 3; i++)
        encode_f64(data, datapos, factor->u.common.z[i]);

    if (factor->u.common.ztruth) {
        encode_u8(data, datapos, 1);

        for (int i = 0; i < 3; i++)
            encode_f64(data, datapos, factor->u.common.ztruth[i]);
    } else {
        encode_u8(data, datapos, 0);
    }

    for (int i = 0; i < 9; i++)
        encode_f64(data, datapos, factor->u.common.W->data[i]);

    april_graph_attr_t *attr = factor->attr;
    stype_encode_object(data, datapos, attr ? attr->stype : NULL, attr);
}

static void *april_graph_factor_xyt_decode(const stype_t *stype, const uint8_t *data, uint32_t *datapos, uint32_t datalen)
{
    int a = decode_u32(data, datapos, datalen);
    int b = decode_u32(data, datapos, datalen);

    double z[3];
    for (int i = 0; i < 3; i++)
        z[i] = decode_f64(data, datapos, datalen);

    double *ztruth = NULL;
    if (decode_u8(data, datapos, datalen)) {
        ztruth = malloc(3*sizeof(double));
        for (int i = 0; i < 3; i++)
            ztruth[i] = decode_f64(data, datapos, datalen);
    }

    matd_t *W = matd_create(3, 3);
    for (int i = 0; i < 9; i++)
        W->data[i] = decode_f64(data, datapos, datalen);

    april_graph_factor_t *factor = april_graph_factor_xyt_create(a, b, z, ztruth, W);
    april_graph_attr_destroy(factor->attr);
    factor->attr = stype_decode_object(data, datapos, datalen, NULL);

    free(ztruth);
    matd_destroy(W);
    return factor;
}

const stype_t stype_april_factor_xyt = { .name = "april_graph_factor_xyt",
                                         .encode = april_graph_factor_xyt_encode,
                                         .decode = april_graph_factor_xyt_decode,
                                         .copy = NULL };

april_graph_factor_t *april_graph_factor_xyt_create(int a, int b, const double *z, const double *ztruth, const matd_t *W)
{
    april_graph_factor_t *factor = calloc(1, sizeof(april_graph_factor_t));

    factor->stype = &stype_april_factor_xyt;
    factor->type = APRIL_GRAPH_FACTOR_XYT_TYPE;
    factor->nnodes = 2;
    factor->nodes = calloc(2, sizeof(int));
    factor->nodes[0] = a;
    factor->nodes[1] = b;
    factor->length = 3;

    factor->copy = xyt_factor_copy;
    factor->eval = xyt_factor_eval;
    factor->state_eval = xyt_factor_state_eval;
    factor->destroy = xyt_factor_destroy;

    factor->u.common.z = doubles_dup(z, 3);
    factor->u.common.ztruth = doubles_dup(ztruth, 3);
    factor->u.common.W = matd_copy(W);

    return factor;
}

/////////////////////////////////////////////////////////////////////////////////////////
// XYT Node
static void xyt_node_update(april_graph_node_t *node, double *dstate)
{
    for (int i = 0; i < 3; i++)
        if(isnan(dstate[i])) return;

    for (int i = 0; i < 3; i++)
        node->state[i] = node->l_point[i] + dstate[i];

    for (int i = 0; i < 3; i++)
        node->delta_X[i] = dstate[i];

    node->state[2] = mod2pi(node->state[2]);
}

static void xyt_node_relinearize(april_graph_node_t *node)
{
    for (int i = 0; i < 3; i++)
        node->l_point[i] = node->state[i];
}

static april_graph_node_t* xyt_node_copy(april_graph_node_t *node)
{
    assert(0);
    return NULL;
}


static void xyt_node_destroy(april_graph_node_t *node)
{
    free(node->state);
    free(node->init);
    free(node->truth);
    free(node->l_point);
    free(node->delta_X);
    if(node->attr) {

        zhash_iterator_t zit;
        april_graph_attr_t *attr = node->attr;
        zhash_iterator_init(attr->hash, &zit);
        char *name;
        struct april_graph_attr_record record;
        while (zhash_iterator_next(&zit, &name, &record)) {
            if (record.stype) {
                record.stype->destroy(record.stype, record.value);
            } else {
                printf("graph attribute %s doesn't have a type\n", name);
            }
        }

        zhash_map_keys(attr->hash, free_key);
        zhash_destroy(attr->hash);
        free(node->attr);
    }
    free(node);
}

static void april_graph_node_xyt_encode(const stype_t *stype, uint8_t *data, uint32_t *datapos, const void *obj)
{
    const april_graph_node_t *node = obj;

    for (int i = 0; i < 3; i++)
        encode_f64(data, datapos, node->state[i]);

    if (node->init) {
        encode_u8(data, datapos, 1);
        for (int i = 0; i < 3; i++)
            encode_f64(data, datapos, node->init[i]);
    } else {
        encode_u8(data, datapos, 0);
    }

    if (node->truth) {
        encode_u8(data, datapos, 1);
        for (int i = 0; i < 3; i++)
            encode_f64(data, datapos, node->truth[i]);
    } else {
        encode_u8(data, datapos, 0);
    }

    april_graph_attr_t *attr = node->attr;
    stype_encode_object(data, datapos, attr ? attr->stype : NULL, attr);
}

static void *april_graph_node_xyt_decode(const stype_t *stype, const uint8_t *data, uint32_t *datapos, uint32_t datalen)
{
    double state[3];

    for (int i = 0; i < 3; i++)
        state[i] = decode_f64(data, datapos, datalen);

    double *init = NULL;
    if (decode_u8(data, datapos, datalen)) {
        init = malloc(3*sizeof(double));
        for (int i = 0; i < 3; i++)
            init[i] = decode_f64(data, datapos, datalen);
    }

    double *truth = NULL;
    if (decode_u8(data, datapos, datalen)) {
        truth = malloc(3*sizeof(double));
        for (int i = 0; i < 3; i++)
            truth[i] = decode_f64(data, datapos, datalen);
    }

    april_graph_node_t *node = april_graph_node_xyt_create(state, init, truth);
    april_graph_attr_destroy(node->attr);
    node->attr = stype_decode_object(data, datapos, datalen, NULL);

    free(init);
    free(truth);
    return node;
}

const stype_t stype_april_node_xyt = { .name = "april_graph_node_xyt",
                                       .encode = april_graph_node_xyt_encode,
                                       .decode = april_graph_node_xyt_decode,
                                       .copy = NULL };

april_graph_node_t *april_graph_node_xyt_create(const double *state, const double *init, const double *truth)
{
    april_graph_node_t *node = calloc(1, sizeof(april_graph_node_t));
    node->type = APRIL_GRAPH_NODE_XYT_TYPE;
    node->length = 3;
    node->state = doubles_dup(state, 3);
    node->init = doubles_dup(init, 3);
    node->truth = doubles_dup(truth, 3);
    node->l_point = doubles_dup(state, 3);
    node->delta_X = calloc(3, sizeof(double));

    node->update = xyt_node_update;
    node->copy = xyt_node_copy;
    node->relinearize = xyt_node_relinearize;
    node->destroy = xyt_node_destroy;

    node->stype = &stype_april_node_xyt;
    return node;
}

void april_graph_xyt_stype_init()
{
    stype_register(&stype_april_factor_xyt);
    stype_register(&stype_april_node_xyt);
}
