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

int main(int argc, char *argv[])
{
    APRILSAM_VERSION();

    setlinebuf(stdout);
    setlinebuf(stderr);

    stype_register_basic_types();
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

        double ninit[3] =  {0, 0, 0};                           // Node initial position
        double nstate[3] = {1, 1, 1 };                          // Node state after optimization
        double ntruth[3] = {2, 2, 2};                           // Ground truth
        april_graph_node_t *node = april_graph_node_xyt_create(nstate,
                                                               ninit,
                                                               ntruth);
        node->UID = zarray_size(graph->nodes);
        zarray_add(graph->nodes, &node);
        printf("============== Save graph ===================\n");
        for(int i=0; i<zarray_size(graph->nodes); i++) {
            april_graph_node_t *node;
            zarray_get(graph->nodes, i, &node);
            printf("node_%d init : (%.2f,%.2f,%.2f) \n", i, node->init[0], node->init[1], node->init[2]);
            printf("node_%d state: (%.2f,%.2f,%.2f) \n", i, node->state[0], node->state[1], node->state[2]);
        }
        april_graph_save(graph, getopt_get_string(gopt, "path"));
        april_graph_destroy(graph);
    }
    {
        //load graph
        april_graph_t *graph = april_graph_create_from_file(getopt_get_string(gopt, "path"));
        printf("============== Load graph ===================\n");
        for(int i=0; i<zarray_size(graph->nodes); i++) {
            april_graph_node_t *node;
            zarray_get(graph->nodes, i, &node);
            printf("node_%d init : (%.2f,%.2f,%.2f) \n", i, node->init[0], node->init[1], node->init[2]);
            printf("node_%d state: (%.2f,%.2f,%.2f) \n", i, node->state[0], node->state[1], node->state[2]);
        }
        april_graph_destroy(graph);
    }

    getopt_destroy(gopt);
}
