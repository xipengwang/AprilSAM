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
   AprilSAM evaluation: ./aprilsam_demo --datapath ../data/M3500.txt
   Batch update only evaluation: ./aprilsam_demo --datapath ../data/M3500.txt --batch_update_only
 */
typedef struct state state_t;
struct state
{
    april_graph_t *graph;
    april_graph_t *loaded_graph;
    struct april_graph_cholesky_param *chol_param;

    double total_time;
    double step_time;

    bool batch_update_only;
    int event_idx; //simulate pose by pose; It is the idx of the pose id in true_path;
};

static int convert_datafile_to_graph(april_graph_t *graph, FILE *f)
{
    char str[20];
    int ID;
    double init[3];
    int aidx;
    int bidx;
    double T[3];
    matd_t *W = matd_create(3, 3);
    while(fscanf(f, "%s", str) != -1) {
        if(!strcmp(str, "VERTEX2")) {
            if(!fscanf(f, "%d %lf %lf %lf", &ID, &init[0], &init[1], &init[2]))
                exit(-1);
            double _state[3] = { init[0], init[1], init[2] };
            double truth[3] = { init[0], init[1], init[2] };
            april_graph_node_t *node = april_graph_node_xyt_create(_state,
                                                                   init,
                                                                   truth);
            zarray_add(graph->nodes, &node);
        } else if(!strcmp(str, "EDGE2")) {
            //IDout   IDin   dx   dy   dth   I11   I12  I22  I33  I13  I23
            if(!fscanf(f, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
                       &aidx, &bidx, &T[0], &T[1], &T[2],
                       &(W->data[0]), &(W->data[1]), &(W->data[4]), &(W->data[8]), &(W->data[2]), &(W->data[5])))
                return -1;

            april_graph_node_t *na, *nb;
            zarray_get(graph->nodes, aidx, &na);
            zarray_get(graph->nodes, bidx, &nb);
            april_graph_factor_t *factor = april_graph_factor_xyt_create(aidx, bidx, T, NULL, W);
            //HACK: iSAM2 w10000 dataset has different format of Manhatan
            if(fabs(bidx - aidx) == 1) {
                april_graph_factor_attr_put(factor, stype_get("string"), "type", strdup("odom"));
            } else {
                april_graph_factor_attr_put(factor, stype_get("string"), "type", strdup("scan"));
            }
            zarray_add(graph->factors, &factor);
        } else {
            printf("error!");
            fclose(f);
            return -1;
        }
    }
    matd_destroy(W);
    printf("%d nodes,  factors: %d \n", zarray_size(graph->nodes), zarray_size(graph->factors));
    fclose(f);
    return 0;
}

static void optimize_cholesky(state_t *state)
{
    int64_t utime0 = utime_now();
    april_graph_cholesky(state->graph, state->chol_param);
    int64_t utime1 = utime_now();
    state->step_time = (utime1 - utime0) / 1.0E3;
    state->total_time += (utime1 - utime0) / 1.0E3;
}

static void optimize_chol_inc(state_t *state)
{
    int64_t utime0 = utime_now();
    april_graph_cholesky_inc(state->graph, state->chol_param);
    int64_t utime1 = utime_now();
    state->step_time = (utime1 - utime0) / 1.0E3;
    state->total_time += (utime1 - utime0) / 1.0E3;
}

static int simulate_on_exist_graph(state_t *state)
{
    printf("Step: %d / %d \n", state->event_idx, zarray_size(state->loaded_graph->nodes));
    if(state->event_idx == zarray_size(state->loaded_graph->nodes))
        return -1;

    april_graph_node_t *loaded_graph_node;
    zarray_get(state->loaded_graph->nodes, state->event_idx , &loaded_graph_node);
    april_graph_node_t *currnode = NULL;
    currnode = april_graph_node_xyt_create(loaded_graph_node->init,
                                           loaded_graph_node->init,
                                           loaded_graph_node->truth);
    zarray_add(state->graph->nodes, &currnode);

    if(state->event_idx == 0) {
        april_graph_factor_t *factor;
        matd_t * W = matd_create_data(3,3, (double []) {10000, 0, 0,
                    0, 10000, 0,
                    0, 0, 1000});
        factor = april_graph_factor_xytpos_create(zarray_size(state->graph->nodes) - 1,
                                                  (double[3]){0,0,0},
                                                  NULL,
                                                  W);
        zarray_add(state->graph->factors, &factor);
        matd_destroy(W);
        state->event_idx++;
        return 0;
    }

    april_graph_factor_t *factor;
    //TODO: could be better
    for (int fidx = 0; fidx < zarray_size(state->loaded_graph->factors); fidx++) {
        april_graph_factor_t *exist_factor;
        zarray_get(state->loaded_graph->factors, fidx, &exist_factor);
        bool add_node_flag = false;
        for (int i = 0; i < exist_factor->nnodes; i++) {
            if(exist_factor->nodes[i] > state->event_idx) {
                add_node_flag = false;
                break;
            }
            if(exist_factor->nodes[i] == state->event_idx) {
                //We may want to add this factor
                add_node_flag = true;
            }
        }
        if(add_node_flag) {
            if(exist_factor->type == APRIL_GRAPH_FACTOR_XYT_TYPE) {
                factor = exist_factor->copy(exist_factor);
                /* printf("Add factor %f,%f,%f\n", */
                /*        factor->u.common.z[0], factor->u.common.z[1], factor->u.common.z[2]); */
                /* matd_print(factor->u.common.W, "%f, "); */
                char *type = april_graph_factor_attr_get(factor, "type");
                if(type) {
                    if(!strcmp(type, "odom")) {
                        //if it is odom factor, update currnode state based on this:
                        double *z = factor->u.common.z;
                        int aidx = factor->nodes[0];
                        int bidx = factor->nodes[1];
                        april_graph_node_t *na, *nb;
                        zarray_get(state->graph->nodes, aidx, &na);
                        zarray_get(state->graph->nodes, bidx, &nb);
                        if(aidx < bidx) {
                            // A + Z = B
                            doubles_xyt_mul(na->state, z, nb->state);
                            nb->relinearize(nb);
                        } else {
                            // A + Z = B
                            // A = B + inv(Z)
                            double inv_z[3];
                            doubles_xyt_inv(z, inv_z);
                            doubles_xyt_mul(nb->state, inv_z, na->state);
                            na->relinearize(na);
                        }
                    }
                    else if(!strcmp(type, "scan")) {
                        //if it is odom factor, update currnode state based on this:
                        //double *z = factor->u.common.z;
                        int aidx = factor->nodes[0];
                        int bidx = factor->nodes[1];
                        //printf("--node id:%d, %d \n", aidx, bidx);
                        april_graph_node_t *na, *nb;
                        zarray_get(state->graph->nodes, aidx, &na);
                        zarray_get(state->graph->nodes, bidx, &nb);
                    }
                    else {
                        int aidx = factor->nodes[0];
                        int bidx = factor->nodes[1];
                        april_graph_node_t *na, *nb;
                        zarray_get(state->graph->nodes, aidx, &na);
                        zarray_get(state->graph->nodes, bidx, &nb);
                    }
                }
                zarray_add(state->graph->factors, &factor);
            }
        }
    }
    state->event_idx++;
    return 0;
}

void simulate_event(state_t *state)
{
    while(1) {
        if(simulate_on_exist_graph(state))
            break;
        if(!state->batch_update_only && state->event_idx > 1) {
            optimize_chol_inc(state);
        } else {
            optimize_cholesky(state);
        }
        double error = april_graph_chi2(state->graph);
        printf("Chi squared error: %f \nStep running time: %.3f ms, Total running time: %.3f ms \n",
                error,
                state->step_time, state->total_time);
    }
}

int main(int argc, char *argv[])
{
    APRILSAM_VERSION();

    setlinebuf(stdout);
    setlinebuf(stderr);

    april_graph_stype_init();
    stype_register_basic_types();
    getopt_t *gopt = getopt_create();
    getopt_add_bool(gopt,   'h',  "help", 0, "Show usage");
    getopt_add_bool(gopt,   '\0', "batch_update_only", 0,  "loaded dataset file path");
    getopt_add_string(gopt, '\0', "datapath",    "",    "loaded dataset file path");
    getopt_add_string(gopt, '\0', "graphpath",  "../data/M3500.graph",   "loaded graph file path");
    getopt_add_int(gopt,    '\0', "nthreshold",  "100", "Batch update if more than nthreshold nodes with significant change");
    getopt_add_double(gopt, '\0', "delta_xy",    "0.1", "re-linearization xy threshold");
    getopt_add_double(gopt, '\0', "delta_theta", "0.1", "re-linearization theta threshold");

    if (!getopt_parse(gopt, argc, argv, 1) || getopt_get_bool(gopt, "help")) {
        getopt_do_usage(gopt);
        return 1;
    }

    state_t *state = calloc(1, sizeof(state_t));
    const char *datapath = getopt_get_string(gopt,"datapath");
    if(!strlen(datapath)) {
        state->loaded_graph = april_graph_create_from_file(getopt_get_string(gopt,"graphpath"));
    }
    else {
        FILE *f = fopen(datapath, "r");
        if(f==NULL) {
            exit(-1);
        }

        state->loaded_graph = april_graph_create();
        if(convert_datafile_to_graph(state->loaded_graph, f)) {
            exit(-1);
        }
        april_graph_save(state->loaded_graph, "/tmp/loaded.graph");
    }

    //Initialize optimizer param
    state->chol_param = calloc(1,sizeof(april_graph_cholesky_param_t));
    april_graph_cholesky_param_init(state->chol_param);

    state->chol_param->show_timing = 0;
    state->chol_param->delta_xy =  getopt_get_double(gopt, "delta_xy");
    state->chol_param->delta_theta = getopt_get_double(gopt, "delta_theta");
    state->chol_param->nthreshold = getopt_get_int(gopt, "nthreshold");
    state->batch_update_only = getopt_get_bool(gopt, "batch_update_only");
    //simulate
    state->graph = april_graph_create();
    simulate_event(state);

    //cleanup
    april_graph_cholesky_param_destory(state->chol_param);
    april_graph_destroy(state->graph);
    getopt_destroy(gopt);

}
