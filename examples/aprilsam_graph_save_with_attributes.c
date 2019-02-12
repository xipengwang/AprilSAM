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

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include "aprilsam/aprilsam.h"
#include "aprilsam/common/getopt.h"

/**
   AprilSAM evaluation: ./aprilsam_graph
 */

// An example of uint64 stype
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

// Example of string stype
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

int main(int argc, char *argv[])
{
    APRILSAM_VERSION();

    setlinebuf(stdout);
    setlinebuf(stderr);

    //Register stype
    april_graph_stype_init();
    stype_register(&stype_uint64);
    stype_register(&stype_string);

    getopt_t *gopt = getopt_create();
    getopt_add_bool(gopt,   'h',  "help", 0, "Show usage");
    getopt_add_string(gopt, '\0', "path", "/tmp/saved.graph", "Graph file path");

    if (!getopt_parse(gopt, argc, argv, 1) || getopt_get_bool(gopt, "help")) {
        getopt_do_usage(gopt);
        return 1;
    }

    {
        //save graph
        april_graph_t *graph = april_graph_create();
        april_graph_attr_put(graph, stype_get("string"), "graph-name", strdup("AprilSAM-Graph"));

        double ninit[3] =  {0, 0, 0};                           // Node initial position
        double nstate[3] = {1, 1, 1 };                          // Node state after optimization
        double ntruth[3] = {2, 2, 2};                           // Ground truth
        //create a xyt node
        april_graph_node_t *node = april_graph_node_xyt_create(nstate,
                                                               ninit,
                                                               ntruth);
        uint64_t *id = calloc(1, sizeof(uint64_t));
        *id = 1024;
        //If user doesn't provide stype, it is user responsibility to free the object(e.g id)
        april_graph_node_attr_put(node, stype_get("uint64"), "ID", id);
        node->UID = zarray_size(graph->nodes);
        zarray_add(graph->nodes, &node);
        //create a xytpos factor
        april_graph_factor_t *factor;
        matd_t * W = matd_create_data(3,3, (double []) {1, 0, 0,
                    0, 1, 0,
                    0, 0, 1});
        factor = april_graph_factor_xytpos_create(0, (double[3]){0,0,0}, NULL, W);
        april_graph_factor_attr_put(factor, stype_get("string"), "factor-type", strdup("geopin"));
        zarray_add(graph->factors, &factor);
        matd_destroy(W);

        printf("============== Save graph ===================\n");
        for(int i=0; i<zarray_size(graph->nodes); i++) {
            april_graph_node_t *node;
            zarray_get(graph->nodes, i, &node);
            printf("Node ID: %ld \n", *id);
            printf("node_%d init : (%.2f,%.2f,%.2f) \n", i, node->init[0], node->init[1], node->init[2]);
            printf("node_%d state: (%.2f,%.2f,%.2f) \n", i, node->state[0], node->state[1], node->state[2]);
        }
        for(int i=0; i<zarray_size(graph->factors); i++) {
            april_graph_factor_t *factor;
            zarray_get(graph->factors, i, &factor);
            printf("factor type: %s \n", "geopin");
        }
        april_graph_save(graph, getopt_get_string(gopt, "path"));
        april_graph_destroy(graph);
    }
    printf("\n");
    {
        //load graph
        april_graph_t *graph = april_graph_create_from_file(getopt_get_string(gopt, "path"));
        char *graph_name = april_graph_attr_get(graph, "graph-name");
        printf("============== Load graph ===================\n");
        printf("Graph name: %s\n", graph_name);
        for(int i=0; i<zarray_size(graph->nodes); i++) {
            april_graph_node_t *node;
            zarray_get(graph->nodes, i, &node);
            uint64_t *id = april_graph_node_attr_get(node, "ID");
            printf("Node ID: %ld \n", *id);
            printf("node_%d init : (%.2f,%.2f,%.2f) \n", i, node->init[0], node->init[1], node->init[2]);
            printf("node_%d state: (%.2f,%.2f,%.2f) \n", i, node->state[0], node->state[1], node->state[2]);
        }
        for(int i=0; i<zarray_size(graph->factors); i++) {
            april_graph_factor_t *factor;
            zarray_get(graph->factors, i, &factor);
            char *factor_type = april_graph_factor_attr_get(factor, "factor-type");
            printf("factor type: %s \n", factor_type);
        }
        april_graph_destroy(graph);
    }

    getopt_destroy(gopt);
}
