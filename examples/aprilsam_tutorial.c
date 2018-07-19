#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include "aprilsam/aprilsam.h"
#include "aprilsam/common/getopt.h"

/**
   AprilSAM evaluation: ./aprilsam_tutorial
 */
typedef struct state state_t;
struct state
{
    april_graph_t *graph;
    struct april_graph_cholesky_param *chol_param;
    double total_time;
    double step_time;
    bool batch_update_only;
};

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

void print_state(state_t *state, int step_id)
{
    printf("\n==================== Step: %d ======================= \n", step_id);
    double error = april_graph_chi2(state->graph);
    printf("Chi squared error: %f \nStep running time: %.3f ms, Total running time: %.3f ms \n",
           error,
           state->step_time, state->total_time);
    //print out state:
    for(int i = 0; i < zarray_size(state->graph->nodes); i++) {
        april_graph_node_t *node;
        zarray_get(state->graph->nodes, i, &node);
        printf("node_%d = {%.2f, %.2f, %.2f} \n", node->UID, node->state[0], node->state[1], node->state[2]);
    }
}

void simulate_event(state_t *state)
{
    //Dog leg dataset with 6 nodes (x, y, heading): (0,0,0), (0,1,0), (0,2,0), (0,3,0), (0,4,0), (0,5,0)
    april_graph_t *graph = state->graph;
    //step 1: Add node_0
    {
        double ninit[3] = {0,0,0};                           // Node initial position
        double nstate[3] = { ninit[0], ninit[1], ninit[2] }; // Node state after optimization
        double ntruth[3] = { ninit[0], ninit[1], ninit[2] }; // Ground truth
        april_graph_node_t *node = april_graph_node_xyt_create(nstate,
                                                               ninit,
                                                               ntruth);
        node->UID = zarray_size(graph->nodes);
        zarray_add(graph->nodes, &node);

        april_graph_factor_t *factor;                                // Add a geopin factor for the first node.
        matd_t * W = matd_create_data(3,3,                           // Information matrix assoicated with the factor
                                      (double []) {10000, 0, 0,
                                              0, 10000, 0,
                                              0, 0, 1000});
        factor = april_graph_factor_xytpos_create(zarray_size(graph->nodes) - 1,
                                                  (double[3]){0,0,0},
                                                  NULL,
                                                  W);
        zarray_add(graph->factors, &factor);
        matd_destroy(W);
        optimize_cholesky(state); //For the first time, you need run a batch update. //TODO: remove this constraint.
        print_state(state, 1);
    }
    //step 2: Add node_1; and add odometry xyt factor between node_0 and node_1
    {
        double ninit[3] = {1,0,0};                           // Node initial position
        double nstate[3] = { ninit[0], ninit[1], ninit[2] }; // Node state after optimization
        double ntruth[3] = { ninit[0], ninit[1], ninit[2] }; // Ground truth
        april_graph_node_t *na;
        april_graph_node_t *nb = april_graph_node_xyt_create(nstate,
                                                             ninit,
                                                             ntruth);
        nb->UID = zarray_size(graph->nodes);
        zarray_get(graph->nodes, zarray_size(graph->nodes)-1, &na);
        zarray_add(graph->nodes, &nb);
        double z[3];
        doubles_xyt_inv_mul(na->init, nb->init, z);
        matd_t *W = matd_create(3,3);
        MATD_EL(W, 0, 0) = 1.0 / pow(0.1, 2);               // Information matrix assoicated with the factor
        MATD_EL(W, 1, 1) = 1.0 / pow(0.1, 2);
        MATD_EL(W, 2, 2) = 1.0 / pow(to_radians(1), 2);
        april_graph_factor_t *factor = april_graph_factor_xyt_create(na->UID, nb->UID, z, NULL, W);
        zarray_add(graph->factors, &factor);
        matd_destroy(W);
        if(state->batch_update_only) {
            optimize_cholesky(state);
        } else {
            optimize_chol_inc(state);
        }
        print_state(state, 2);
    }
    //step 3: Add node_2; and add odometry xyt factor between node_1 and node_2
    {
        double ninit[3] = {2,0,0};                           // Node initial position
        double nstate[3] = { ninit[0], ninit[1], ninit[2] }; // Node state after optimization
        double ntruth[3] = { ninit[0], ninit[1], ninit[2] }; // Ground truth
        april_graph_node_t *na;
        april_graph_node_t *nb = april_graph_node_xyt_create(nstate,
                                                             ninit,
                                                             ntruth);
        zarray_get(graph->nodes, zarray_size(graph->nodes)-1, &na);
        nb->UID = zarray_size(graph->nodes);
        zarray_add(graph->nodes, &nb);
        double z[3];
        doubles_xyt_inv_mul(na->init, nb->init, z);
        matd_t *W = matd_create(3,3);
        MATD_EL(W, 0, 0) = 1.0 / pow(0.1, 2);
        MATD_EL(W, 1, 1) = 1.0 / pow(0.1, 2);
        MATD_EL(W, 2, 2) = 1.0 / pow(to_radians(1), 2);
        april_graph_factor_t *factor = april_graph_factor_xyt_create(na->UID, nb->UID, z, NULL, W);
        zarray_add(graph->factors, &factor);
        matd_destroy(W);
        if(state->batch_update_only) {
            optimize_cholesky(state);
        } else {
            optimize_chol_inc(state);
        }
        print_state(state, 3);
    }
    //step 4: Add node_3; and add odometry xyt factor between node_2 and node_3
    {
        double ninit[3] = {3,0,0};                           // Node initial position
        double nstate[3] = { ninit[0], ninit[1], ninit[2] }; // Node state after optimization
        double ntruth[3] = { ninit[0], ninit[1], ninit[2] }; // Ground truth
        april_graph_node_t *na;
        april_graph_node_t *nb = april_graph_node_xyt_create(nstate,
                                                             ninit,
                                                             ntruth);
        zarray_get(graph->nodes, zarray_size(graph->nodes)-1, &na);
        nb->UID = zarray_size(graph->nodes);
        zarray_add(graph->nodes, &nb);
        double z[3];
        doubles_xyt_inv_mul(na->init, nb->init, z);
        matd_t *W = matd_create(3,3);
        MATD_EL(W, 0, 0) = 1.0 / pow(0.1, 2);
        MATD_EL(W, 1, 1) = 1.0 / pow(0.1, 2);
        MATD_EL(W, 2, 2) = 1.0 / pow(to_radians(1), 2);
        april_graph_factor_t *factor = april_graph_factor_xyt_create(na->UID, nb->UID, z, NULL, W);
        zarray_add(graph->factors, &factor);
        matd_destroy(W);
        if(state->batch_update_only) {
            optimize_cholesky(state);
        } else {
            optimize_chol_inc(state);
        }
        print_state(state, 4);
    }
    //step 5: Add node_4; and add odometry xyt factor between node_3 and node_4
    {
        double ninit[3] = {4,0,0};                           // Node initial position
        double nstate[3] = { ninit[0], ninit[1], ninit[2] }; // Node state after optimization
        double ntruth[3] = { ninit[0], ninit[1], ninit[2] }; // Ground truth
        april_graph_node_t *na;
        april_graph_node_t *nb = april_graph_node_xyt_create(nstate,
                                                             ninit,
                                                             ntruth);
        zarray_get(graph->nodes, zarray_size(graph->nodes)-1, &na);
        nb->UID = zarray_size(graph->nodes);
        zarray_add(graph->nodes, &nb);
        double z[3];
        doubles_xyt_inv_mul(na->init, nb->init, z);
        matd_t *W = matd_create(3,3);
        MATD_EL(W, 0, 0) = 1.0 / pow(0.1, 2);
        MATD_EL(W, 1, 1) = 1.0 / pow(0.1, 2);
        MATD_EL(W, 2, 2) = 1.0 / pow(to_radians(1), 2);
        april_graph_factor_t *factor = april_graph_factor_xyt_create(na->UID, nb->UID, z, NULL, W);
        zarray_add(graph->factors, &factor);
        matd_destroy(W);
        if(state->batch_update_only) {
            optimize_cholesky(state);
        } else {
            optimize_chol_inc(state);
        }
        print_state(state, 5);
    }
    //step 6.a: Add node_5; and add odometry xyt factor between node_4 and node_5
    {
        double ninit[3] = {5,0,0};                           // Node initial position
        double nstate[3] = { ninit[0], ninit[1], ninit[2] }; // Node state after optimization
        double ntruth[3] = { ninit[0], ninit[1], ninit[2] }; // Ground truth
        april_graph_node_t *na;
        april_graph_node_t *nb = april_graph_node_xyt_create(nstate,
                                                             ninit,
                                                             ntruth);
        zarray_get(graph->nodes, zarray_size(graph->nodes)-1, &na);
        nb->UID = zarray_size(graph->nodes);
        zarray_add(graph->nodes, &nb);
        double z[3];
        doubles_xyt_inv_mul(na->init, nb->init, z);
        matd_t *W = matd_create(3,3);
        MATD_EL(W, 0, 0) = 1.0 / pow(0.1, 2);
        MATD_EL(W, 1, 1) = 1.0 / pow(0.1, 2);
        MATD_EL(W, 2, 2) = 1.0 / pow(to_radians(1), 2);
        april_graph_factor_t *factor = april_graph_factor_xyt_create(na->UID, nb->UID, z, NULL, W);
        zarray_add(graph->factors, &factor);
        matd_destroy(W);
    }
    //step 6.b Add a xyt factor between node_0 and node_5. Node_0 thinks node_5 is at (5,1,0).
    {
        april_graph_node_t *na;
        april_graph_node_t *nb;
        zarray_get(graph->nodes, 0, &na);
        zarray_get(graph->nodes, zarray_size(graph->nodes)-1, &nb);
        double z[3];
        double tmp[3] = {5, 1, 0};     // The first node think the last node is at (5,1,0) not (5,0,0);
        doubles_xyt_inv_mul(na->init, tmp, z);
        matd_t *W = matd_create(3,3);
        MATD_EL(W, 0, 0) = 1.0 / pow(0.1, 2);
        MATD_EL(W, 1, 1) = 1.0 / pow(0.1, 2);
        MATD_EL(W, 2, 2) = 1.0 / pow(to_radians(1), 2);
        april_graph_factor_t *factor = april_graph_factor_xyt_create(na->UID, nb->UID, z, NULL, W);
        zarray_add(graph->factors, &factor);
        matd_destroy(W);
        if(state->batch_update_only) {
            optimize_cholesky(state);
        } else {
            optimize_chol_inc(state);
        }
        print_state(state, 6);
    }
}

int main(int argc, char *argv[])
{
    APRILSAM_VERSION();

    setlinebuf(stdout);
    setlinebuf(stderr);

    stype_register_basic_types();
    april_graph_stype_init();
    getopt_t *gopt = getopt_create();
    getopt_add_bool(gopt,   'h',  "help", 0, "Show usage");
    getopt_add_int(gopt,    '\0', "nthreshold",  "100", "Batch update if more than nthreshold nodes with significant change");
    getopt_add_double(gopt, '\0', "delta_xy",    "0.1", "re-linearization xy threshold");
    getopt_add_double(gopt, '\0', "delta_theta", "0.1", "re-linearization theta threshold");
    getopt_add_bool(gopt,   '\0', "batch_update_only", 0,  "loaded dataset file path");

    if (!getopt_parse(gopt, argc, argv, 1) || getopt_get_bool(gopt, "help")) {
        getopt_do_usage(gopt);
        return 1;
    }

    state_t *state = calloc(1, sizeof(state_t));

    //Initialize optimizer param
    state->chol_param = calloc(1,sizeof(april_graph_cholesky_param_t));
    april_graph_cholesky_param_init(state->chol_param);
    //Don't show timeprofiling of aprilsam
    state->chol_param->show_timing = 0;
    //re-linearization xy threshold(meter)
    state->chol_param->delta_xy =  getopt_get_double(gopt, "delta_xy");
    //re-linearization theta threshold(radian)
    state->chol_param->delta_theta = getopt_get_double(gopt, "delta_theta");
    //Batch update if more than nthreshold nodes with significant at least delta_xy or delta_theta changes.
    state->chol_param->nthreshold = getopt_get_int(gopt, "nthreshold");
    state->batch_update_only = getopt_get_bool(gopt, "batch_update_only");

    //simulate
    state->graph = april_graph_create();
    simulate_event(state);

    //cleanup
    april_graph_cholesky_param_destory(state->chol_param);
    april_graph_destroy(state->graph);
    getopt_destroy(gopt);
    free(state);
}
