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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <limits.h>

#include "aprilsam.h"

void APRILSAM_VERSION()
{
    printf("================\n");
    printf("| APRILSAM 0.0 |\n");
    printf("================\n");
    printf("\n");
}

/** POSE Graph for xyt
   convert matrix column id to graph node id Assume all node just xyt!!!!! */
static int *heap_minimum_degree_ordering(smatd_t *mat, april_graph_t *g);

void april_graph_cholesky_param_init(april_graph_cholesky_param_t *param)
{

    memset(param, 0, sizeof(april_graph_cholesky_param_t));
    param->ordering = NULL;
    param->tikhanov = 0.0001;
    param->show_timing = 0;
    param->nreordering = 1;
    param->chol = NULL;
    param->A = NULL;
    param->tr = NULL;
    if(param->delta_x) {
        free(param->delta_x);
        param->delta_x = NULL;
    }
    if(param->ordering) {
        free(param->ordering);
        param->ordering = NULL;
    }
}

void april_graph_cholesky_param_destory(april_graph_cholesky_param_t *param)
{
    if(!param)
        return;
    if(param->chol)
        smatd_chol_destroy(param->chol);
    if(param->A)
        smatd_destroy(param->A);
    if(param->tr)
        search_tree_destroy(param->tr);
    if(param->delta_x)
        free(param->delta_x);
    if(param->B)
        free(param->B);
    if(param->y)
        free(param->y);
    if(param->ordering)
        free(param->ordering);
    free(param);
}

void april_graph_cholesky(april_graph_t *graph, april_graph_cholesky_param_t *_param)
{
    // nothing to do
    if (zarray_size(graph->nodes) == 0 || zarray_size(graph->factors) == 0)
        return;

    if (_param->nreordering) {

        april_graph_cholesky_param_t *param = _param;

        timeprofile_t *tp = timeprofile_create();
        timeprofile_stamp(tp, "begin");

        int *allocated_ordering = NULL; // we'll free this one.
        int *use_ordering = param->ordering; // this could be from .param or one we create.

        // make symbolic matrix for variable nreordering.
        smatd_t *Asym = smatd_create(zarray_size(graph->nodes), zarray_size(graph->nodes));
        for (int fidx = 0; fidx < zarray_size(graph->factors); fidx++) {
            april_graph_factor_t *factor;
            zarray_get(graph->factors, fidx, &factor);

            for (int i = 0; i < factor->nnodes; i++) {
                for (int j = 0; j < factor->nnodes; j++) {
                    smatd_set(Asym, factor->nodes[i], factor->nodes[j], 1);
                }
            }
        }

        timeprofile_stamp(tp, "make symbolic");

        //check filling
        if(1) {
            //allocated_ordering = exact_minimum_degree_ordering(Asym);
            allocated_ordering = heap_minimum_degree_ordering(Asym, graph);
        } else {
            allocated_ordering = calloc(zarray_size(graph->nodes), sizeof(int));
            for (int i = 0; i < zarray_size(graph->nodes); i++) {
                allocated_ordering[i] = i;
            }
        }
        timeprofile_stamp(tp, "compute ordering");

        //relinearize the point
        for (int i = 0; i < zarray_size(graph->nodes); i++) {
            april_graph_node_t *node;
            zarray_get(graph->nodes, i, &node);
            node->relinearize(node);
        }

        use_ordering = allocated_ordering;
        smatd_destroy(Asym);

        // idxs[j]: what index in x do the state variables for node j start at?
        int *idxs = calloc(zarray_size(graph->nodes), sizeof(int));
        int xlen = 0;
        for (int i = 0; i < zarray_size(graph->nodes); i++) {
            april_graph_node_t *node;
            zarray_get(graph->nodes, use_ordering[i], &node);
            idxs[use_ordering[i]] = xlen;
            xlen += node->length;
        }
        assert(xlen > 0);

        // we'll solve normal equations, Ax = B
        smatd_t *A = smatd_create(xlen, xlen);
        double  *B = calloc(xlen, sizeof(double));
        for (int i = 0; i < zarray_size(graph->factors); i++) {
            april_graph_factor_t *factor;
            zarray_get(graph->factors, i, &factor);
            april_graph_factor_eval_t *eval = factor->eval(factor, graph, NULL);

            for (int z0 = 0; z0 < factor->nnodes; z0++) {
                int n0 = factor->nodes[z0];

                matd_t *JatW = matd_op("M'*M", eval->jacobians[z0], eval->W);

                for (int z1 = 0; z1 < factor->nnodes; z1++) {
                    int n1 = factor->nodes[z1];

                    matd_t *JatWJb = matd_op("M*M", JatW, eval->jacobians[z1]);

                    for (int row = 0; row < JatWJb->nrows; row++) {
                        for (int col = 0; col < JatWJb->ncols; col++) {
                            if(row+idxs[n0] > col+idxs[n1])
                                continue;
                            int a = row+idxs[n0];
                            int b = col+idxs[n1];
                            double v = smatd_get(A, a, b);
                            double new_v = v + MATD_EL(JatWJb, row, col);
                            smatd_set(A, a, b, new_v);
                        }
                    }

                    matd_destroy(JatWJb);
                }

                matd_t *R = matd_create_data(eval->length, 1, eval->r);
                matd_t *JatWr = matd_op("M*M", JatW, R);
                for (int row = 0; row < JatWr->nrows; row++)
                    B[idxs[n0]+row] += MATD_EL(JatWr, row, 0);

                matd_destroy(R);
                matd_destroy(JatW);
                matd_destroy(JatWr);
            }

            april_graph_factor_eval_destroy(eval);
        }

        if (param->tikhanov > 0) {
            //TODO: A is upper triangle, we can A->rows[i].values[0]
            double lambda = param->tikhanov;
            for (int i = 0; i < xlen; i++) {
                double v = smatd_get(A, i, i);
                smatd_set(A, i, i, v + lambda);
            }
        }
        timeprofile_stamp(tp, "build A, B");

        //HACK
        /* xlen = 3; */
        /* A = smatd_create(xlen, xlen); */
        /* smatd_set(A, 0, 0, 6);smatd_set(A, 0, 1, 15);smatd_set(A, 0, 2, 55); */
        /* smatd_set(A, 1, 1, 55);smatd_set(A, 1, 2, 225); */
        /* smatd_set(A, 2, 2, 979); */

        //TODO: A is upper triangle, we probably don't need to use smatd_upper_right(A)
        //Record chol result
        cs *CS_A = cs_spalloc (xlen, xlen, 1, 1, 1);
        for(int k = 0; k < xlen; k++) {
            svecd_t *v = &(A->rows[k]);
            for(int j = 0; j < v->nz; j++) {
                int col = v->indices[j];
                cs_entry(CS_A, k, col, v->values[j]);
                if(k!=col)
                    cs_entry(CS_A, col, k, v->values[j]);
            }
        }

        smatd_t *u = smatd_create(xlen, xlen);
        {
            cs *A = cs_triplet(CS_A);
            int order = -1;
            css *S ;
            csn *N ;
            S = cs_schol (A, order) ;		/* ordering and symbolic analysis */
            N = cs_chol (A, S) ;		/* numeric Cholesky factorization */
            //convert N to chol->u
            cs *L = N->L;
            for(int k = 0; k < xlen; k++) {
                svecd_t *v = &(u->rows[k]);
                //L is column compressed
                int row_b = L->p[k];
                int row_s = L->p[k+1];
                int mincap = row_s - row_b;
                v->alloc = mincap;
                v->indices = calloc(1,mincap*sizeof(int));
                v->values = calloc(1,mincap*sizeof(double));
                v->nz = mincap;
                memcpy(v->indices, &(L->i[row_b]), mincap*sizeof(int));
                memcpy(v->values, &(L->x[row_b]), mincap*sizeof(double));
            }
            //chol->u = u;
            cs_sfree (S) ;
            cs_nfree (N) ;
            cs_spfree(A);
        }
        //smatd_destroy(u);
        smatd_chol_t *chol = calloc(1, sizeof(smatd_chol_t));
        chol->u = u;
        timeprofile_stamp(tp, "convert chol");
        //chol = smatd_chol(A);
        if(param->chol) {
            smatd_chol_destroy(param->chol);
        }
        param->chol = chol;
        timeprofile_stamp(tp, "chol decompose");
        if(param->tr) {
            search_tree_destroy(param->tr);
            param->tr = NULL;
        }
        param->tr = search_tree_create_from_smat(chol->u, allocated_ordering, idxs, graph);
        param->tr->delta_xy = param->delta_xy;
        param->tr->delta_theta = param->delta_theta;
        timeprofile_stamp(tp, "generate tree");

        if(param->B) {
            free(param->B);
        }
        param->B = B;

        if(param->A) {
            smatd_destroy(param->A);
        }
        param->A = A;
        if(param->ordering)
            free(param->ordering);

        param->ordering = allocated_ordering;
        param->nreordering = zarray_size(graph->nodes);
        param->factor_num = zarray_size(graph->factors);
        {
            timeprofile_t *tp = timeprofile_create();
            timeprofile_stamp(tp, "begin");

            if(param->y) {
                free(param->y);
            }
            param->y = calloc(param->chol->u->ncols, sizeof(double));
            double *x = calloc(param->chol->u->ncols, sizeof(double));
            smatd_chol_solve_full(param->chol, param->B, param->y, x);
            int *use_ordering = param->ordering;
            int *idxs = calloc(param->nreordering, sizeof(int));
            int xlen = 0;

            assert(param->nreordering <= zarray_size(graph->nodes));
            for (int i = 0; i < param->nreordering; i++) {
                april_graph_node_t *node;
                zarray_get(graph->nodes, use_ordering[i], &node);
                idxs[use_ordering[i]] = xlen;
                xlen += node->length;
            }

            for (int i = param->nreordering - 1; i >= 0; i--) {
                april_graph_node_t *node;
                zarray_get(graph->nodes, i, &node);
                node->update(node, &x[idxs[i]]);
            }
            timeprofile_stamp(tp, "solve");
            if (param->show_timing)
                timeprofile_display(tp);
            timeprofile_destroy(tp);
            //CSPARSE
            if(0){
                timeprofile_t *tp = timeprofile_create();
                timeprofile_stamp(tp, "begin---CSPARSE");
                cs *A = cs_triplet(CS_A);
                timeprofile_stamp(tp, "Change A,B");
                double *b = B;
                int order = 0;
                double *x ;
                css *S ;
                csn *N ;
                int n, ok ;
                n = A->n ;
                S = cs_schol (A, order) ;		/* ordering and symbolic analysis */
                N = cs_chol (A, S) ;		/* numeric Cholesky factorization */
                x = cs_malloc (n, sizeof (double)) ;
                timeprofile_stamp(tp, "chol");
                ok = (S && N && x) ;
                if (ok)
                {
                    cs_ipvec (n, S->Pinv, b, x) ;	/* x = P*b */
                    cs_lsolve (N->L, x) ;		/* x = L\x */
                    cs_ltsolve (N->L, x) ;		/* x = L'\x */
                    cs_pvec (n, S->Pinv, x, b) ;	/* b = P'*x */
                }
                cs_free (x) ;
                cs_sfree (S) ;
                cs_nfree (N) ;
                for (int i = param->nreordering - 1; i >= 0; i--) {
                    april_graph_node_t *node;
                    zarray_get(graph->nodes, i, &node);
                    node->update(node, &B[idxs[i]]);
                }
                cs_spfree(A);
                free(B);
                timeprofile_stamp(tp, "solve-2");
                if (param->show_timing)
                    timeprofile_display(tp);
                timeprofile_destroy(tp);
            }
            free(x);
            free(idxs);
            cs_spfree(CS_A);
            if(param->delta_x) {
                free(param->delta_x);
                param->delta_x = x;
            }

        }
        free(idxs);
        timeprofile_destroy(tp);
    }
    else {
        assert(0);
    }
}

void april_graph_cholesky_inc(april_graph_t *graph, april_graph_cholesky_param_t *param)
{
    // nothing to do
    if (zarray_size(graph->nodes) == 0 || zarray_size(graph->factors) == 0)
        return;
    if(!param->chol)
        return;
    if(param->factor_num == zarray_size(graph->factors))
        return;

    // Use chol instead of long variable name param->chol
    smatd_chol_t *chol = param->chol;
    timeprofile_t *tp = timeprofile_create();
    timeprofile_stamp(tp, "begin");

    // Ordering is already augmented, same size as number of nodes
    int *use_ordering = realloc(param->ordering, sizeof(int)*zarray_size(graph->nodes));
    for(int i = param->nreordering; i < zarray_size(graph->nodes); i++) {
        use_ordering[i] = i;
    }
    param->ordering = use_ordering;

    // idxs[j]: what index in x do the state variables for node j start at?
    int *idxs = calloc(zarray_size(graph->nodes), sizeof(int));
    int xlen = 0;
    for (int i = 0; i < zarray_size(graph->nodes); i++) {
        april_graph_node_t *node;
        zarray_get(graph->nodes, use_ordering[i], &node);
        idxs[use_ordering[i]] = xlen;
        xlen += node->length;
    }
    assert(xlen > 0);
    //timeprofile_stamp(tp, "augment ordering");
    // Augment reward B, y, and chol->u if number of nodes gets larger
    if(zarray_size(graph->nodes) > param->nreordering) {
        //TODO: Could be smart if we pre-allocate lots of memory of A, B, etc...
        double *B = calloc(xlen, sizeof(double));
        memcpy(B, param->B, sizeof(double)*chol->u->nrows);
        free(param->B);
        param->B = B;

        double *y = calloc(xlen, sizeof(double));
        memcpy(y, param->y, sizeof(double)*chol->u->nrows);
        free(param->y);
        param->y= y;

        svecd_t *rows = calloc(xlen, sizeof(svecd_t));
        memcpy(rows, chol->u->rows, sizeof(svecd_t) * chol->u->nrows);
        free(chol->u->rows);
        chol->u->rows = rows;

        //XXX: Why use realloc not working?
        //chol->u->rows = realloc(chol->u->rows, xlen * sizeof(svecd_t));
        //param->B = realloc(param->B, sizeof(double) * xlen);

        chol->u->nrows = xlen;
        chol->u->ncols = xlen;

        for (int i = 0; i < chol->u->nrows; i++) {
            chol->u->rows[i].length = chol->u->ncols;
        }
        //XXX: This part is used for keeping previous information matrix; it seems not helping a lot for current
        //testing, and will investigate this later.
        svecd_t *A_rows = calloc(xlen, sizeof(svecd_t));
        memcpy(A_rows, param->A->rows, sizeof(svecd_t) * param->A->nrows);
        free(param->A->rows);
        param->A->rows = A_rows;
        param->A->nrows = xlen;
        param->A->ncols = xlen;

        for (int i = 0; i < chol->u->nrows; i++) {
            param->A->rows[i].length = param->A->ncols;
        }
    }

    // Augment the search tree, but keep the root node unchanged since we haven't updated new nodes.
    search_tree_t *tr = param->tr;
    assert(param->tr);
    int nnodes = tr->nnodes;
    tr->nnodes = zarray_size(graph->nodes);
    int root_id = tr->root->id;
    april_graph_node_t **nodes_array = (april_graph_node_t**)graph->nodes->data;
    if (tr->nnodes >= tr->nalloc) {
        tr->nalloc = tr->nnodes * 2;
        search_tree_node_t *tmp = calloc(tr->nalloc, sizeof(search_tree_node_t));
        memcpy(tmp, tr->nodes, nnodes * sizeof(search_tree_node_t));
        free(tr->nodes);
        tr->nodes = tmp;
        tr->linearized_nodes = realloc(tr->linearized_nodes, tr->nalloc * sizeof(int));
    }
    tr->root = &(tr->nodes[use_ordering[root_id]]);
    for(int i = nnodes; i < tr->nnodes; i++) {
        tr->nodes[i].nalloc = 8;
        tr->nodes[i].nchildren = 0;
        tr->nodes[i].children = calloc(tr->nodes[i].nalloc, sizeof(int));
        tr->nodes[i].parent = -1;
        tr->nodes[i].g_node = nodes_array[i];
        tr->nodes[i].g_node->UID = i;
        tr->nodes[i].id = i;
        assert(tr->root->id == root_id);
    }

    // Figure out affected nodes and mark their 'label_changed' = 1, we assume node are xyt_node.
    // Tree has x, y, t as nodes, so  id % 3 == 0 is x, id % 3 == 1 is y, id % 3 == 2 is t.
    // x's parent is always y of same graph node, y's parent is always t, and t's parent is another graph node's x.
    tr->naffected = 0;
    for (int i = param->factor_num; i < zarray_size(graph->factors); i++) {
        april_graph_factor_t *factor;
        zarray_get(graph->factors, i, &factor);
        for (int z0 = 0; z0 < factor->nnodes; z0++) {
            int n0 = factor->nodes[z0];
            search_tree_node_t *node= &(param->tr->nodes[n0]);
            while(!node->label_changed) {
                node->label_changed = 1;
                tr->naffected++;
                if(node->parent != -1)
                    node = &(param->tr->nodes[node->parent]);
                else
                    break;
            }
        }
    }
    //timeprofile_stamp(tp, "check modified nodes");

    //timeprofile_stamp(tp, "augment data");
    smatd_chol_reconstruct_tr(chol, param->ordering, param->tr, param->y);
    timeprofile_stamp(tp, "Reconstruct");

    // we'll solve normal equations, Ax = B
    double  *B = param->B;
    //Add new information into A and B.
    for (int i = param->factor_num; i < zarray_size(graph->factors); i++) {
        april_graph_factor_t *factor;
        zarray_get(graph->factors, i, &factor);
        april_graph_factor_eval_t *eval = factor->eval(factor, graph, NULL);
        for (int z0 = 0; z0 < factor->nnodes; z0++) {
            int n0 = factor->nodes[z0];
            matd_t *JatW = matd_op("M'*M", eval->jacobians[z0], eval->W);
            for (int z1 = 0; z1 < factor->nnodes; z1++) {
                int n1 = factor->nodes[z1];
                matd_t *JatWJb = matd_op("M*M", JatW, eval->jacobians[z1]);
                for (int row = 0; row < JatWJb->nrows; row++) {
                    for (int col = 0; col < JatWJb->ncols; col++) {
                        if(row+idxs[n0] > col+idxs[n1])
                            continue;
                        double v = smatd_get(chol->u, row+idxs[n0], col+idxs[n1]);
                        smatd_set(chol->u, row+idxs[n0], col+idxs[n1], v + MATD_EL(JatWJb, row, col));
                        v = smatd_get(param->A, row+idxs[n0], col+idxs[n1]);
                        smatd_set(param->A, row+idxs[n0], col+idxs[n1], v + MATD_EL(JatWJb, row, col));
                    }
                }
                matd_destroy(JatWJb);
            }
            matd_t *R = matd_create_data(eval->length, 1, eval->r);
            matd_t *JatWr = matd_op("M*M", JatW, R);
            for (int row = 0; row < JatWr->nrows; row++) {
                B[idxs[n0]+row] += MATD_EL(JatWr, row, 0);
                param->y[idxs[n0]+row] += MATD_EL(JatWr, row, 0);
                assert(idxs[n0]+row < xlen);
            }
            matd_destroy(R);
            matd_destroy(JatW);
            matd_destroy(JatWr);
        }
        april_graph_factor_eval_destroy(eval);
    }
    // We have param->factor_num factors
    param->factor_num = zarray_size(graph->factors);
    timeprofile_stamp(tp, "build A, B");
    // Incremental cholesky decomposition
    // TODO: should we do partial reordering here before we perform cholesky decomposition
    smatd_chol_inc_tr(chol, param->tr);
    timeprofile_stamp(tp, "chol decompose");
    search_tree_append(param->tr, chol->u, param->ordering, idxs);
    timeprofile_stamp(tp, "append tree");

    if (param->show_timing) {
        timeprofile_display(tp);
    }
    //HACK: takes longer than 150 ms
    if(1.0E-3 * timeprofile_total_utime(tp) > param->batch_time / 3) {
        param->tr->start_over = INT_MAX;
    }
    timeprofile_destroy(tp);

    param->nreordering = zarray_size(graph->nodes);
    april_graph_cholesky_inc_solver(graph, param, idxs);
    free(idxs);

    if(param->tr->start_over > param->nthreshold) {
        free(param->ordering);
        param->ordering = NULL;
        int64_t utime0 = utime_now();
        april_graph_cholesky(graph, param);
        int64_t utime1 = utime_now();
        param->batch_time = (utime1 - utime0) / 1.0E3;
        param->tr->start_over = 0;
        param->tr->nlinearized_nodes = 0;
    }
}

void april_graph_cholesky_inc_solver(april_graph_t *graph, april_graph_cholesky_param_t *param, int *idxs)
{
    if(param->nreordering) {
        timeprofile_t *tp = timeprofile_create();
        timeprofile_stamp(tp, "begin");
        double *x = calloc(param->chol->u->ncols, sizeof(double));
        smatd_chol_solve_tr(param->chol, param->B, param->y, x, param->tr);
        assert(param->nreordering <= zarray_size(graph->nodes));
        timeprofile_stamp(tp, "solve");
        if (param->show_timing)
            timeprofile_display(tp);
        timeprofile_destroy(tp);
        if(param->delta_x) {
            free(param->delta_x);
            param->delta_x = x;
        } else {
            free(x);
        }
    }
}

search_tree_t *search_tree_create(int nnodes)
{
    search_tree_t *tr = calloc(1, sizeof(search_tree_t));
    tr->nnodes = nnodes;
    tr->nalloc = tr->nnodes * 2;
    tr->nodes = calloc(tr->nalloc, sizeof(search_tree_node_t));
    for(int i = 0; i < tr->nalloc; i++) {
        tr->nodes[i].nalloc = 8;
        tr->nodes[i].children = calloc(tr->nodes[i].nalloc, sizeof(int));
        tr->nodes[i].parent = -1;
    }
    return tr;
}

search_tree_t *search_tree_create_from_smat(smatd_t *u, int *ordering, int *idxs, april_graph_t *g)
{
    april_graph_node_t **nodes_array = (april_graph_node_t**)g->nodes->data;
    search_tree_t *tr = calloc(1, sizeof(search_tree_t));
    tr->nnodes = u->nrows / 3;
    tr->nalloc = tr->nnodes;
    tr->nodes = calloc(tr->nalloc, sizeof(search_tree_node_t));
    tr->nlinearized_nodes = 0;
    tr->linearized_nodes = calloc(tr->nalloc, sizeof(int));
    for(int i = 0; i < tr->nnodes; i++) {
        tr->nodes[i].nalloc = 8;
        tr->nodes[i].nchildren = 0;
        tr->nodes[i].children = calloc(tr->nodes[i].nalloc, sizeof(int));
        tr->nodes[i].parent = -1;
        tr->nodes[i].g_node = nodes_array[i];
        tr->nodes[i].g_node->UID = i;
    }
    /* for(int i = tr->nnodes; i < tr->nalloc; i++) { */
    /*     tr->nodes[i].nalloc = 8; */
    /*     tr->nodes[i].children = calloc(tr->nodes[i].nalloc, sizeof(int)); */
    /*     tr->nodes[i].parent = -1; */
    /* } */
    for(int ui = tr->nnodes-2; ui >= 0; ui--) {
        int child_node_id = ordering[ui];
        search_tree_node_t *child_node = &(tr->nodes[child_node_id]);
        svecd_t *urowi = &u->rows[idxs[child_node_id] + 2]; // Hard Code as xyt_node, so 2 here.
        if(urowi->nz > 1) {
            int parent_node_id = ordering[urowi->indices[1] / 3]; // Hard Code as xyt_node, so 3 here.
            search_tree_node_t *parent_node = &(tr->nodes[parent_node_id]);
            child_node->id = ui;
            child_node->parent = parent_node_id;
            if(parent_node->nchildren == parent_node->nalloc) {
                parent_node->nalloc = parent_node->nalloc * 2;
                parent_node->children = realloc(parent_node->children,
                                                parent_node->nalloc * sizeof(int));
            }
            parent_node->children[parent_node->nchildren] = child_node_id;
            parent_node->nchildren++;
        }
    }
    int ui = tr->nnodes-1;
    tr->root = &(tr->nodes[ordering[ui]]);
    tr->root->id = ui;
    return tr;
}

static void print_node(search_tree_t *tr, search_tree_node_t* node)
{
    printf("node id: %d, nchildren:%d, label:%d,%d ",
           node->g_node->UID, node->nchildren, node->label_changed, node->label_relinearized);

    if(node->parent != -1) {
        printf("parent %d ->", tr->nodes[node->parent].g_node->UID);
    }

    for (int i = 0; i < node->nchildren; i++) {
        printf("%d,", tr->nodes[node->children[i]].g_node->UID);
    }
    printf("\n");
    for (int i = 0; i < node->nchildren; i++) {
        print_node(tr, &(tr->nodes[node->children[i]]));
    }
}

void search_tree_print(search_tree_t *tr)
{
    printf("root id: %d, nnodes:%d\n", tr->root->g_node->UID, tr->nnodes);
    /* for(int i = 0; i < tr->nnodes; i++) { */
    /*     search_tree_node_t *node = &(tr->nodes[i]); */
    /*     printf(":%d=>", node->id); */
    /*     for (int i = 0; i < node->nchildren; i++) { */
    /*         printf("%d,", node->children[i]->id); */
    /*     } */
    /*     printf("\n"); */
    /* } */
    /* return; */
    print_node(tr, tr->root);
}

void search_tree_destroy(search_tree_t *tr)
{
    for(int i = 0; i < tr->nalloc; i++) {
        free(tr->nodes[i].children);
    }
    free(tr->nodes);
    free(tr->linearized_nodes);
    free(tr);
}

static void recursive_smatd_ltransposetriangle_solve_tr(search_tree_t *tr, search_tree_node_t *node,
                                                        smatd_t *u, double *x)
{
    if(!node->label_changed)
        return;
    int k = node->id * 3; //HARD CODE: xyt_node, so 3 here.
    for (int i = 0; i < node->nchildren; i++) {
        recursive_smatd_ltransposetriangle_solve_tr(tr, &(tr->nodes[node->children[i]]), u, x);
    }
    for (int i = k; i <= k+2; i++) {  //HARD CODE: xyt_node, so 2 here.
        svecd_t *urowi = &u->rows[i];
        x[i] = x[i] / urowi->values[0];
        for (int pos = 1; pos < urowi->nz; pos++) {
            int uidx = urowi->indices[pos];
            x[uidx] -= x[i] * urowi->values[pos];
        }
    }
}

static void solve_node(search_tree_t *tr, search_tree_node_t* node, smatd_t *u, double *b, double *x)
{
    int k = node->id * 3;  //HARD CODE: xyt_node, so 3 here.
    for (int i = 2+k; i >= k; i--) { //HARD CODE: xyt_node, so 2 here.
        double bi = b[i];
        svecd_t *urowi = &u->rows[i];
        // this row should end with a value along the diagonal.
        // NB: rank-deficient matrices could have a zero here?
        assert(urowi->nz > 0);
        assert(urowi->indices[0] == i);

        TYPE diag = urowi->values[0];

        // compute dot product of urowi and the last (i-1) elements of
        // y.  Note: exploit fact that u is exactly triangular-- the
        // first element will be the diagonal (which we want to skip)
        for (int pos = 1; pos < urowi->nz ; pos++)
            bi -= urowi->values[pos]*x[urowi->indices[pos]];

        x[i] = bi / diag;
        if(i % 3 == 0) {
            if(fabs(x[i]) > tr->delta_xy || fabs(x[i+1]) > tr->delta_xy || fabs(x[i+2]) > tr->delta_theta) {
                if(!node->label_relinearized) {
                    node->label_relinearized = 1;
                    tr->linearized_nodes[tr->nlinearized_nodes++] = node->g_node->UID;
                    tr->start_over += 1;
                    tr->total_delta_xy += (fabs(x[i]) + fabs(x[i+1]));
                } else {
                    tr->total_delta_xy += (fabs(x[i]) + fabs(x[i+1]) - node->g_node->delta_X[0] - node->g_node->delta_X[1]);
                }
            }
            node->g_node->delta_X[0] = x[i];
            node->g_node->delta_X[1] = x[i+1];
            node->g_node->delta_X[2] = x[i+2];
            if(tr->naffected > 5) {
                node->label_changed = 0;
            } else {
                if(node->label_changed == 1) {
                    node->label_changed = 0;
                } else {
                    //check if x = x[i], y = x[i+1], theta = x[i+2] only have small change with respect to previous changes
                    //if so, just return;
                    double delta_delta_x = node->g_node->delta_X[0] - x[i];
                    double delta_delta_y = node->g_node->delta_X[1] - x[i+1];
                    double delta_delta_theta = node->g_node->delta_X[2] - x[i+2];
                    //HARD CODE THRESHOLD
                    if(fabs(delta_delta_x) < 0.01 &&
                       fabs(delta_delta_y) < 0.01 && fabs(delta_delta_theta) < 2 * M_PI / 180) {
                        return;
                    }
                }
            }
            node->g_node->update(node->g_node, &x[i]);
        }
    }
    for (int i = 0; i < node->nchildren; i++) {
        solve_node(tr, &(tr->nodes[node->children[i]]), u, b, x);
    }
}

// Solve the system Ax=b. User provides storage for x.
void smatd_chol_solve_tr(smatd_chol_t *chol, const double *b, double *y, double *x, search_tree_t *tr)
{
    // U'Ux = b
    // U'y = b
    smatd_t *u =chol->u;
    recursive_smatd_ltransposetriangle_solve_tr(tr, tr->root, u, y);
    solve_node(tr, tr->root, chol->u, y, x);
}

static void recursive_reconstruct_y(smatd_chol_t *chol, double *y, search_tree_t *tr, search_tree_node_t *node)
{
    if(!node->label_changed)
        return;
    smatd_t *u = chol->u;
    int k = node->id * 3;  //HARD CODE: xyt_node, so 3 here.
    for (int i = 2+k; i >= k; i--) { //HARD CODE: xyt_node, so 2 here.
        //reconstruct y
        svecd_t *urowi = &u->rows[i];
        for (int pos = 1; pos < urowi->nz; pos++) {
            int uidx = urowi->indices[pos];
            y[uidx] += y[i] * urowi->values[pos];
        }
        y[i] = y[i] * urowi->values[0];
    }
    for (int i = 0; i < node->nchildren; i++) {
        recursive_reconstruct_y(chol, y, tr, &(tr->nodes[node->children[i]]));
    }
}

static void recursive_reconstruct_u(smatd_chol_t *chol, search_tree_t *tr, search_tree_node_t *node)
{
    if(!node->label_changed)
        return;
    smatd_t *u = chol->u;
    int k = node->id * 3;  //HARD CODE: xyt_node, so 3 here.
    for (int i = 2+k; i >= k; i--) {  //HARD CODE: xyt_node, so 2 here.
        //reconstruct chol-u
        svecd_t *urowi = &u->rows[i];
        double d = svecd_get(urowi, i);
        for (int pos = 0; pos < urowi->nz; pos++) {
            int j = urowi->indices[pos];
            if (j > i) {
                // apply an update to a row below us
                double s = urowi->values[pos];
                svecd_t *urowj = &u->rows[j];
                svecd_add_i0_x(urowj, urowj, urowi, s, 0, pos);
            }
        }
        svecd_scale(urowi, d);
    }
    for (int i = 0; i < node->nchildren; i++) {
        recursive_reconstruct_u(chol, tr, &(tr->nodes[node->children[i]]));
    }
}

void smatd_chol_reconstruct_tr(smatd_chol_t *chol, int *ordering, search_tree_t *tr, double *y)
{
    //call reconstruct_y before reconstruct_u
    timeprofile_t *tp = timeprofile_create();
    //timeprofile_stamp(tp, "begin");
    recursive_reconstruct_y(chol, y, tr, tr->root);
    //timeprofile_stamp(tp, "y");
    recursive_reconstruct_u(chol, tr, tr->root);
    //timeprofile_stamp(tp, "u");
    timeprofile_display(tp);
    timeprofile_destroy(tp);
}

static void recursive_chol_decomp(smatd_chol_t *chol, search_tree_t *tr, search_tree_node_t *node)
{
    if(!node->label_changed)
        return;
    for (int i = 0; i < node->nchildren; i++) {
        recursive_chol_decomp(chol, tr, &(tr->nodes[node->children[i]]));
    }
    int k = node->id * 3;  //HARD CODE: xyt_node, so 3 here.
    for (int i = k; i <= k+2; i++) {  //HARD CODE: xyt_node, so 2 here.
        smatd_t *u = chol->u;
        svecd_t *urowi = &u->rows[i];
        double d = svecd_get(urowi, i);
        d = sqrt(d);

        svecd_scale(urowi, 1.0 / d);
        for (int pos = 0; pos < urowi->nz; pos++) {
            int j = urowi->indices[pos];
            if (j > i) {
                // apply an update to a row below us
                double s = urowi->values[pos];
                svecd_t *urowj = &u->rows[j];
                svecd_add_i0_x(urowj, urowj, urowi, -s, 0, pos);
            }
        }
    }
}

void smatd_chol_inc_tr(smatd_chol_t *chol, search_tree_t *tr)
{
    //chol->u should be upper right already.
    /* smatd_t *u = smatd_upper_right(chol->u); */
    /* smatd_destroy(chol->u); */
    /* chol->u = u */
    smatd_t *u = chol->u;
    recursive_chol_decomp(chol, tr, tr->root);
    int n = u->nrows / 3;  //HARD CODE: xyt_node, so 3 here.
    int is_spd = 1;
    for (int k = tr->root->id+1; k < n; k++) {
        for (int i = 3*k; i <= 3*k+2; i++) {  //HARD CODE: xyt_node, so 2 here.
            svecd_t *urowi = &u->rows[i];
            double d = svecd_get(urowi, i);
            is_spd &= (d > 0);
            d = sqrt(d);
            svecd_scale(urowi, 1.0 / d);
            for (int pos = 0; pos < urowi->nz; pos++) {
                int j = urowi->indices[pos];
                if (j > i) {
                    // apply an update to a row below us
                    double s = urowi->values[pos];
                    svecd_t *urowj = &u->rows[j];
                    svecd_add_i0_x(urowj, urowj, urowi, -s, 0, pos);
                }
            }
        }
    }
    chol->is_spd = is_spd;
}

static void recursive_tree_append(search_tree_t *tr, search_tree_node_t *node, smatd_t *u, int *ordering, int *idxs)
{
    if(!node->label_changed)
        return;
    int ui = node->id;
    int child_node_id = ordering[ui];
    search_tree_node_t *child_node = &(tr->nodes[child_node_id]);
    svecd_t *urowi = &u->rows[idxs[child_node_id] + 2]; // Hard Code as xyt_node, so 2 here.
    if(urowi->nz > 1) {
        int parent_node_id = ordering[urowi->indices[1] / 3]; // Hard Code as xyt_node, so 3 here.
        search_tree_node_t *parent_node = &(tr->nodes[parent_node_id]);
        child_node->id = ui;
        //We don't want to add duplicated child_node;
        if(child_node->parent != -1) {
            if(child_node->parent == parent_node_id) {
                goto next_node;
            } else {
                //TODO:Bugs in here?
                //parent change, we need to remove this child from its previous parent.
                search_tree_node_t *old_parent_node = &(tr->nodes[child_node->parent]);
                int label_id = old_parent_node->nchildren;
                for(int i=0; i<old_parent_node->nchildren; i++) {
                    search_tree_node_t *c_node = &(tr->nodes[old_parent_node->children[i]]);
                    if(c_node->id == child_node->id) {
                        assert(c_node->g_node->UID == child_node->g_node->UID);
                        label_id = i;
                        break;
                    }
                }
                old_parent_node->nchildren--;
                for(int i = label_id+1; i<old_parent_node->nchildren; i++) {
                    old_parent_node->children[i-1] = old_parent_node->children[i];
                }
            }
        }
        child_node->parent = parent_node_id;
        if(parent_node->nchildren >= parent_node->nalloc) {
            parent_node->nalloc = parent_node->nchildren * 2;
            parent_node->children = realloc(parent_node->children, parent_node->nalloc*sizeof(int));
        }
        parent_node->children[parent_node->nchildren] = child_node_id;
        parent_node->nchildren++;
    }

  next_node:
    for (int i = 0; i < node->nchildren; i++) {
        recursive_tree_append(tr, &(tr->nodes[node->children[i]]), u, ordering, idxs);
    }
}

void search_tree_append(search_tree_t *tr, smatd_t *u, int *ordering, int *idxs)
{
    recursive_tree_append(tr, tr->root, u, ordering, idxs);
    int root_id = tr->root->id;
    for (int ui = u->nrows / 3 - 2; ui > root_id; ui--) {
        int child_node_id = ordering[ui];
        search_tree_node_t *child_node = &(tr->nodes[child_node_id]);
        svecd_t *urowi = &u->rows[idxs[child_node_id] + 2]; // Hard Code as xyt_node, so 2 here.
        if(urowi->nz > 1) {
            int parent_node_id = ordering[urowi->indices[1] / 3]; // Hard Code as xyt_node, so 3 here.
            search_tree_node_t *parent_node = &(tr->nodes[parent_node_id]);
            child_node->id = ui;
            if(child_node->parent != -1) {
                if(child_node->parent == parent_node_id) {
                    continue;
                }
            }
            child_node->parent = parent_node_id;
            if(parent_node->nchildren >= parent_node->nalloc) {
                parent_node->nalloc = parent_node->nchildren * 2;
                parent_node->children = realloc(parent_node->children, parent_node->nalloc*sizeof(int));
            }
            parent_node->children[parent_node->nchildren] = child_node_id;
            parent_node->nchildren++;
        }
    }
    int ui = u->nrows / 3 - 1;
    tr->root = &(tr->nodes[ordering[ui]]);
    tr->root->id = ui;
}

struct node
{
    int *neighbors;
    int nneighbors;
    int alloc;
    int id;
    int degree;
};
typedef struct node node_t;

static int *heap_minimum_degree_ordering(smatd_t *mat, april_graph_t *g)
{
    //We could use a array instead of hashing.
    zhash_t *nodes_hash = zhash_create(sizeof(uint32_t), sizeof(zqueue_t*), zhash_uint32_hash, zhash_uint32_equals);
    zmaxheap_t *nodes_heap = zmaxheap_create(sizeof(zqueue_t*));

    int nnodes = mat->nrows;
    struct node *nodes = calloc(nnodes, sizeof(struct node));

    for (int rowi = 0; rowi < mat->nrows; rowi++) {
        svecd_t *vec = &mat->rows[rowi];
        nodes[rowi].alloc = vec->nz * 2 + 8;
        nodes[rowi].neighbors = malloc(nodes[rowi].alloc * sizeof(int));
        nodes[rowi].id = rowi;
        nodes[rowi].nneighbors = 0;
        for (int i = 0; i < vec->nz; i++) {
            if (vec->indices[i] == rowi)
                continue;
            nodes[rowi].neighbors[nodes[rowi].nneighbors++] = vec->indices[i];
        }
    }
    int *set_marker = calloc(nnodes, sizeof(int));
    {
        //most recent node set the lowest score.
        int rowi =  mat->nrows - 1;
        uint32_t key = nodes[rowi].nneighbors + 2 * rowi;
        zqueue_t *node_list;
        if(zhash_get(nodes_hash, &key, &node_list)) {
            node_t *node_p = &(nodes[rowi]);
            zqueue_add_back(node_list, &node_p);
        } else {
            node_list = zqueue_create(sizeof(node_t*));
            node_t *node_p = &(nodes[rowi]);
            zqueue_add_back(node_list, &node_p);
            zhash_put(nodes_hash, &key, &node_list, NULL, NULL);
            zmaxheap_add(nodes_heap, &node_list, -1.0 * key);
        }
        set_marker[rowi] = 1;
        /* for(int i = 0; i < nodes[rowi].nneighbors; i++) { */
        /*     int choose_idx = nodes[rowi].neighbors[i]; */
        /*     for (int idx = choose_idx - 5; idx < choose_idx + 5; idx++) { */
        /*         if(idx < 0 || idx > mat->nrows-1) */
        /*             continue; */
        /*         if(set_marker[idx] == 1) */
        /*             continue; */
        /*         uint32_t key = nodes[idx].nneighbors + rowi; */
        /*         zqueue_t *node_list; */
        /*         if(zhash_get(nodes_hash, &key, &node_list)) { */
        /*             node_t *node_p = &(nodes[idx]); */
        /*             zqueue_add_back(node_list, &node_p); */
        /*         } else { */
        /*             node_list = zqueue_create(sizeof(node_t*)); */
        /*             node_t *node_p = &(nodes[idx]); */
        /*             zqueue_add_back(node_list, &node_p); */
        /*             zhash_put(nodes_hash, &key, &node_list, NULL, NULL); */
        /*             zmaxheap_add(nodes_heap, &node_list, -1.0 * key); */
        /*         } */
        /*         set_marker[idx] = 1; */
        /*     } */
        /* } */

        for(int i = 0; i < nodes[rowi].nneighbors; i++) {
            int choose_idx = nodes[rowi].neighbors[i];
            for (int idx = choose_idx - 5; idx < choose_idx + 5; idx++) {
                if(idx < 0 || idx > mat->nrows-1)
                    continue;
                if(set_marker[idx] == 1)
                    continue;
                uint32_t key = nodes[idx].nneighbors + rowi;
                zqueue_t *node_list;
                if(zhash_get(nodes_hash, &key, &node_list)) {
                    node_t *node_p = &(nodes[idx]);
                    zqueue_add_back(node_list, &node_p);
                } else {
                    node_list = zqueue_create(sizeof(node_t*));
                    node_t *node_p = &(nodes[idx]);
                    zqueue_add_back(node_list, &node_p);
                    zhash_put(nodes_hash, &key, &node_list, NULL, NULL);
                    zmaxheap_add(nodes_heap, &node_list, -1.0 * key);
                }
                set_marker[idx] = 1;
                for(int j = 0; j < nodes[idx].nneighbors; j++) {
                    if(set_marker[j] == 1)
                        continue;
                    uint32_t key = nodes[j].nneighbors + rowi;
                    zqueue_t *node_list;
                    if(zhash_get(nodes_hash, &key, &node_list)) {
                        node_t *node_p = &(nodes[j]);
                        zqueue_add_back(node_list, &node_p);
                    } else {
                        node_list = zqueue_create(sizeof(node_t*));
                        node_t *node_p = &(nodes[j]);
                        zqueue_add_back(node_list, &node_p);
                        zhash_put(nodes_hash, &key, &node_list, NULL, NULL);
                        zmaxheap_add(nodes_heap, &node_list, -1.0 * key);
                    }
                }
            }
        }
    }

    for (int rowi = 0; rowi < mat->nrows - 1; rowi++) {
        if(set_marker[rowi] == 1)
            continue;

        uint32_t key = nodes[rowi].nneighbors;
        zqueue_t *node_list;
        if(zhash_get(nodes_hash, &key, &node_list)) {
            node_t *node_p = &(nodes[rowi]);
            zqueue_add_back(node_list, &node_p);
        } else {
            node_list = zqueue_create(sizeof(node_t*));
            node_t *node_p = &(nodes[rowi]);
            zqueue_add_back(node_list, &node_p);
            zhash_put(nodes_hash, &key, &node_list, NULL, NULL);
            zmaxheap_add(nodes_heap, &node_list, -1.0 * key);
        }
    }

    free(set_marker);
    int orderingi = 0;
    int *ordering = calloc(nnodes, sizeof(int));
    // a hash set used to detect duplicate nodes when marginalizing
    // out nodes. See below.
    int *set = calloc(nnodes, sizeof(int));
    int *marker = calloc(nnodes, sizeof(int));
    int settoken = 0;
    zqueue_t *node_list;
    float v;
    while(zmaxheap_remove_max(nodes_heap, &node_list, &v)) {
        //printf("list size: %d \n" , zqueue_size(node_list));
        //printf("value : %f \n", v);
        while(zqueue_size(node_list)) {
            node_t *node_p;
            zqueue_pop(node_list, &node_p);
            if(marker[node_p->id] ==1) {
                //already eliminated before
                continue;
            }
            //TODO: AMD? calculate degree d instead of using nneighbors?
            if(node_p->nneighbors <= -v) {
                //remove the node
                ordering[orderingi] = node_p->id;
                marker[node_p->id] = 1;
                //printf("id, key, nneighbors: %d, %f, %d \n" , node_p->id, v, node_p->nneighbors);
                orderingi++;
                int bestnodei = node_p->id;
                //updates its neighbors;
                // marginalize-out bestnodei: every neighbor of bestnodei becomes connected.
                for (int aidx = 0; aidx < nodes[bestnodei].nneighbors; aidx++) {
                    int naidx = nodes[bestnodei].neighbors[aidx];

                    struct node *na = &nodes[naidx];

                    // na is the "destination" to which we will add all of the
                    // neighbors of bestnodei.

                    // make the 'set' contain the current neighbors.
                    // set membership is denoted by set[i] == settoken.
                    // we get "free" set tables by just incrementing settoken.
                    settoken++;

                    for (int i = 0; i < na->nneighbors; i++) {

                        // while we're iterating, delete the bestnodei neighbor
                        if (na->neighbors[i] == bestnodei) {
                            // shuffle delete
                            na->neighbors[i] = na->neighbors[na->nneighbors - 1];
                            na->nneighbors--;
                            i--;
                            continue;
                        }

                        set[na->neighbors[i]] = settoken;
                    }

                    // setting these set elements at this point will prevent us
                    // from adding these nodes below.
                    set[bestnodei] = settoken;
                    set[naidx] = settoken;

                    // now loop over the neighbors to add.
                    for (int bidx = 0; bidx < nodes[bestnodei].nneighbors; bidx++) {

                        int nbidx = nodes[bestnodei].neighbors[bidx];

                        // set already contains this member
                        if (set[nbidx] == settoken)
                            continue;

                        // need to add nj to neighbors of ni.
                        if (na->nneighbors + 1 >= na->alloc) {
                            // need to grow the neighbors array.
                            na->alloc *= 2;
                            assert(na->nneighbors > 0 && na->alloc > 0);
                            na->neighbors = realloc(na->neighbors, na->alloc*sizeof(int));
                        }

                        // add the neighbor
                        na->neighbors[na->nneighbors++] = nbidx;
                    }
                    //TODO: Use degree instead of use nneighbors since we modify
                    //degree based on distance score.
                    /* if(na->nneighbors <= -v) { */
                    /*     //if less than key, we add back to this queue, but we need to remember to remove this na from */
                    /*     //higher degree queue by using a marker */
                    /*     zqueue_add_back(node_list, &na); */
                    /*     continue; */
                    /* } */

                    /* if(na->nneighbors < na->degree) { */
                    /*     uint32_t key = na->nneighbors; */
                    /*     na->degree = na->nneighbors; */
                    /*     assert(key != 0); */
                    /*     zqueue_t *node_list; */
                    /*     if(zhash_get(nodes_hash, &key, &node_list)) { */
                    /*         zqueue_add_back(node_list, &na); */
                    /*     } else { */
                    /*         node_list = zqueue_create(sizeof(node_t*)); */
                    /*         zqueue_add_back(node_list, &na); */
                    /*         zhash_put(nodes_hash, &key, &node_list, NULL, NULL); */
                    /*         zmaxheap_add(nodes_heap, &node_list, -1.0 * key); */
                    /*     } */
                    /* } */
                }
            } else {
                //re-insert the node to correct queue;
                uint32_t key = node_p->nneighbors;
                zqueue_t *node_list;
                if(zhash_get(nodes_hash, &key, &node_list)) {
                    zqueue_add_back(node_list, &node_p);
                } else {
                    node_list = zqueue_create(sizeof(node_t*));
                    zqueue_add_back(node_list, &node_p);
                    zmaxheap_add(nodes_heap, &node_list, -1.0 * key);
                }
            }
        }
        zqueue_destroy(node_list);
    }

    for (int rowi = 0; rowi < mat->nrows; rowi++)
        free(nodes[rowi].neighbors);

    free(marker);
    free(nodes);
    free(set);
    zhash_destroy(nodes_hash);
    zmaxheap_destroy(nodes_heap);
    return ordering;
}
