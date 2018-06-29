/* $COPYRIGHT_UM
   $LICENSE_LGPL
*/

#ifndef _APRIL_GRAPH_H
#define _APRIL_GRAPH_H

#include <stdlib.h>

#include "./common/doubles.h"
#include "./common/matd.h"
#include "./common/math_util.h"
#include "./common/smatd.h"
#include "./common/string_util.h"
#include "./common/stype.h"
#include "./common/timeprofile.h"
#include "./common/zarray.h"
#include "./common/zhash.h"
#include "./common/zqueue.h"
#include "./common/zmaxheap.h"


#include "./csparse/csparse.h"

void APRILSAM_VERSION();

/** april_graph_attr */
struct april_graph_attr_record
{
    const stype_t *stype;
    void          *value;
};

typedef struct april_graph_attr april_graph_attr_t;
struct april_graph_attr
{
    zhash_t *hash;
    const stype_t *stype;
};

april_graph_attr_t *april_graph_attr_create();
void april_graph_attr_destroy(april_graph_attr_t *attr);

/** april_graph */
typedef struct april_graph april_graph_t;
struct april_graph
{
    zarray_t *factors;
    zarray_t *nodes;

    april_graph_attr_t *attr;
    const stype_t *stype;
};

/** april_graph_factor_eval */
typedef struct april_graph_factor_eval april_graph_factor_eval_t;
struct april_graph_factor_eval
{
    double chi2;

    // One jacobian for every node that this factor is connected to.
    // (order given by the factor's 'nodes' array). Should be NULL
    // terminated so that this structure can be deallocated without knowing
    // which factor created it.
    matd_t **jacobians;

    int length;
    double *r; // residual (length x 1)
    matd_t *W; // information matrix of observation (length x length)
};

#define APRIL_GRAPH_FACTOR_XYT_TYPE 1
#define APRIL_GRAPH_FACTOR_XYTPOS_TYPE 2

#define APRIL_GRAPH_NODE_XYT_TYPE 100

/** april_graph_factor */

typedef struct april_graph_factor april_graph_factor_t;
struct april_graph_factor
{
    int type; // a unique identifier

    int nnodes;
    int *nodes;

    int length; // how many DOF?

    april_graph_attr_t *attr;

    april_graph_factor_t* (*copy)(april_graph_factor_t *factor);

    // evaluate this factor potential given the current graph. The
    // "eval" object can be expensive to build; for efficiency, a
    // caller can recycle an eval object previously constructed by
    // this object. If non-null, this object MUST be recycled and the
    // factor will return the passed-in "eval" object (with recomputed
    // numerical values). Otherwise, a new "eval" object is created
    // and returned.
    april_graph_factor_eval_t* (*eval)(april_graph_factor_t *factor, april_graph_t *graph, april_graph_factor_eval_t *eval);
    april_graph_factor_eval_t* (*state_eval)(april_graph_factor_t *factor, april_graph_t *graph, april_graph_factor_eval_t *eval);

    void (*destroy)(april_graph_factor_t *factor);

    // storage for implementations. They can use them however they
    // like.
    union {
        struct {
            double *z;
            double *ztruth;
            matd_t *W;
            void   *impl;
        } common;

        struct {
            april_graph_factor_t **factors;
            double *logw;
            int nfactors;
        } max;

        struct {
            void *impl;
        } impl;
    } u;

    const stype_t *stype;
};

/////////////////////////////////////////////////////////////
// april_graph_node

typedef struct april_graph_node april_graph_node_t;
struct april_graph_node
{
    int UID;
    int type;  // a unique identifier
    int length; // # of DOF

    double *state; // length x 1  vectors
    double *init;
    double *truth;
    double *l_point; //linearization point
    double *delta_X; //last delta change

    april_graph_attr_t *attr;  // string (char*) to arbitrary pointer

    april_graph_node_t* (*copy)(april_graph_node_t *node);

    // called after some optimization; generically, should
    // implement state += dstate.
    void (*update)(april_graph_node_t *node, double *dstate);

    void (*relinearize)(april_graph_node_t *node);

    void (*destroy)(april_graph_node_t *node);

    void *impl;

    const stype_t *stype;
};

/////////////////////////////////////////////////////////////
// API

april_graph_t *april_graph_create();
april_graph_t *april_graph_create_from_file(const char *path);
void april_graph_destroy(april_graph_t *graph);

void april_graph_factor_eval_destroy(april_graph_factor_eval_t *eval);

typedef struct search_tree_node search_tree_node_t;
struct search_tree_node
{
    int *children;
    int parent;
    int nalloc;
    int nchildren;
    int id;            // sorted column id
    april_graph_node_t *g_node;
    int label_changed; // if 1, the nodes are affected by adding new factors.
    int label_relinearized;
};

typedef struct search_tree search_tree_t;
struct search_tree
{
    int nnodes;
    int nalloc;
    search_tree_node_t *root;
    search_tree_node_t *nodes;
    int start_over;
    int nlinearized_nodes;
    int *linearized_nodes;
    int isam1_cnt;
    int naffected;
    double delta_xy;
    double delta_theta;
    double total_delta_xy;
    double total_delta_theta;
};
search_tree_t *search_tree_create(int nnodes);
void search_tree_print(search_tree_t *tr);
void search_tree_destroy(search_tree_t *tr);

search_tree_t *search_tree_create_from_smat(smatd_t *u, int *ordering, int *idxs, april_graph_t *g);
void search_tree_append(search_tree_t *tr, smatd_t *u, int *reordering, int *idxs);
void smatd_chol_solve_tr(smatd_chol_t *chol, const double *b, double *y, double *x, search_tree_t *tr);
void smatd_chol_inc_tr(smatd_chol_t *chol, search_tree_t *tr);
void smatd_chol_reconstruct_tr(smatd_chol_t *chol, int *ordering, search_tree_t *tr, double *y);

//Incremental Cholesky
typedef struct april_graph_cholesky_param april_graph_cholesky_param_t;
struct april_graph_cholesky_param
{
    // if non-zero, add tikahnov*I to the information matrix.
    double tikhanov;

    smatd_chol_t *chol; // Record last iteration cholesky result for incremental smoothing

    int factor_num;     // Record number of factors in last iteration.

    // Use the specified node ordering to reduce fill-in. If not
    // specified, an ordering is computed automatically.
    int *ordering;       // Record ordering in last iteration
    int nreordering;     // non-zero enables incremental smoothing, 0 use default cholesky
    int show_timing;     // Boolean; non-zero enable showing time profiling.

    //R^T * R * delta_x = B
    //R^T * y = B
    //R * delta_x = y
    double *delta_x;     // Record delta_x values in last iteration
    double *B;           // Record reward in last iteration
    double *y;           // Record last intermediate reward;
    smatd_t *A;

    search_tree_t *tr;

    double l_thresh;    //check if need to be relinearized
    double delta_thresh;
    int nthreshold; //Batch update if more than nthreshold nodes with significant change

    double batch_time;

    double delta_xy;
    double delta_theta;
};

// initialize to default values.
void april_graph_cholesky_param_init(april_graph_cholesky_param_t *param);
void april_graph_cholesky_param_destory(april_graph_cholesky_param_t *param);
// Compute a Gauss-Newton update on the graph, using the specified
// node ordering. The ordering should specify the order for each
// april_graph_node_t; if NULL, a default ordering is computed. The
// ordering passed in belongs to the caller.
void april_graph_cholesky(april_graph_t *graph, april_graph_cholesky_param_t *param);
void april_graph_cholesky_inc(april_graph_t *graph, april_graph_cholesky_param_t *param);
void april_graph_cholesky_inc_solver(april_graph_t *graph, april_graph_cholesky_param_t *param, int *idxs);
void april_graph_cholesky_qr_inc(april_graph_t *graph, april_graph_cholesky_param_t *param);


int april_graph_dof(april_graph_t *graph);
double april_graph_chi2(april_graph_t *graph);

april_graph_node_t *april_graph_node_xyt_create(const double *state, const double *init, const double *truth);

april_graph_factor_t *april_graph_factor_xyt_create(int a, int b, const double *z, const double *ztruth, const matd_t *W);
april_graph_factor_t *april_graph_factor_xytpos_create(int a, double *z, double *ztruth, matd_t *W);

void april_graph_attr_put(april_graph_t *graph, const stype_t *type, const char *key, void *data);
void* april_graph_attr_get(april_graph_t *graph, const char * key);

void april_graph_factor_attr_put(april_graph_factor_t *factor, const stype_t *type, const char *key, void *data);
void* april_graph_factor_attr_get(april_graph_factor_t *factor, const char * key);

void april_graph_node_attr_put(april_graph_node_t *node, const stype_t *type, const char *key, void *data);
void* april_graph_node_attr_get(april_graph_node_t *node, const char * key);

int april_graph_save(april_graph_t *graph, const char *path);

void april_graph_stype_init();

#endif
