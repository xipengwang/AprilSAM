/* $COPYRIGHT_UM
$LICENSE_LGPL
*/

/*$LICENSE*/

// sparse mat
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>

#include "assert.h"
#include "timeprofile.h"

#include "smatd.h"

// make future smatf/smati implementations easier.
#define TYPE double

#define SVEC_MIN_CAPACITY 16

/////////////////////////////////////////////////
// vec operations

svecd_t *svecd_create(int length)
{
    svecd_t *v = calloc(1, sizeof(svecd_t));
    v->length = length;
    v->indices = NULL;
    v->values = NULL;
    return v;
}

// beware: smatd allocates svecd_ts directly. Do not call
void svecd_destroy(svecd_t *v)
{
    if (!v)
        return;
    if(!v->indices)
        free(v->indices);
    if(!v->values)
        free(v->values);

    // this is why you can't call svecd_destroy on a svec allocated by a smatd
    free(v);
}

void svecd_ensure_capacity(svecd_t *v, int mincap)
{
    if (v->alloc >= mincap)
        return;

    int newcap = v->alloc;
    if (newcap < SVEC_MIN_CAPACITY)
        newcap = SVEC_MIN_CAPACITY;

    while (newcap < mincap) {
        newcap *= 2;
    }

    v->indices = realloc(v->indices, sizeof(int)*newcap);
    v->values = realloc(v->values, sizeof(TYPE)*newcap);
    v->alloc = newcap;
}

// make the ith element correspond to the (idx, v) tuple, moving
// anything after it as necessary.l Position 'i' must be the correct
// position for idx. 'idx' must not already be in the vec.
void svecd_insert(svecd_t *v, int i, int idx, TYPE val)
{
    if (val==0)
        return;

    svecd_ensure_capacity(v, v->nz + 1);

    memmove(&v->indices[i+1], &v->indices[i], sizeof(int)*(v->nz - i));
    memmove(&v->values[i+1], &v->values[i], sizeof(TYPE)*(v->nz - i));

    v->indices[i] = idx;
    v->values[i] = val;
    v->nz++;
}

void svecd_set(svecd_t *v, int idx, TYPE val)
{
    if (v->setcursor < 0 || v->setcursor >= v->nz)
        v->setcursor = v->nz / 2;

    if (v->nz == 0) {
        if (val == 0)
            return;

        svecd_ensure_capacity(v, 1);

        v->indices[0] = idx;
        v->values[0] = val;
        v->nz = 1;
        return;
    }

    if (v->indices[v->setcursor] == idx) {
        v->values[v->setcursor] = val;
        return;
    }

    // search.
    if (v->indices[v->setcursor] < idx) {
        // search up
        while (v->setcursor+1 < v->nz && v->indices[v->setcursor+1] <= idx)
            v->setcursor++;

        if (v->indices[v->setcursor] == idx) {
            v->values[v->setcursor] = val;
            return;
        }

        svecd_insert(v, v->setcursor+1, idx, val);
        return;
    } else {
        // search down
        while (v->setcursor-1 >= 0 && v->indices[v->setcursor-1] >= idx)
            v->setcursor--;

        if (v->indices[v->setcursor] == idx) {
            v->values[v->setcursor] = val;
            return;
        }

        svecd_insert(v, v->setcursor, idx, val);
        return;
    }
}

TYPE svecd_get(svecd_t *v, int idx)
{
    if (v->nz == 0)
        return 0;

    if (v->getcursor >= v->nz || v->getcursor < 0)
        v->getcursor = v->nz / 2;

    if (v->indices[v->getcursor] < idx) {
        // search up
        while (v->getcursor+1 < v->nz && v->indices[v->getcursor+1] <= idx)
            v->getcursor++;

        if (v->indices[v->getcursor]==idx)
            return v->values[v->getcursor];

        return 0;
    } else {
        // search down
        while (v->getcursor-1 >= 0 && v->indices[v->getcursor-1] >= idx)
            v->getcursor--;

        if (v->indices[v->getcursor] == idx)
            return v->values[v->getcursor];

        return 0;
    }

    assert(0);
}

TYPE svecd_dot_product(svecd_t *a, svecd_t *b)
{
    if (a->nz == 0 || b->nz == 0)
        return 0;

    int aidx = 0, bidx = 0;
    int ai = a->indices[aidx], bi = b->indices[bidx];

    TYPE acc = 0;

    while (aidx < a->nz && bidx < b->nz) {

        if (ai == bi) {
            acc += a->values[aidx]*b->values[bidx];
            aidx++;
            bidx++;
            ai = a->indices[aidx];
            bi = b->indices[bidx];
            continue;
        }

        if (ai < bi) {
            aidx++;
            ai = a->indices[aidx];
        } else {
            bidx++;
            bi = b->indices[bidx];
        }
    }

    return acc;
}

// v = v * s
void svecd_scale(svecd_t *v, double scale)
{
    for (int i = 0; i < v->nz; i++) {
        v->values[i] *= scale;
    }
}

// beginning at element i0, scale the element
void svecd_scale_i0(svecd_t *v, double scale, int i0)
{
    for (int i = 0; i < v->nz; i++) {
        if (v->indices[i] >= i0)
            v->values[i] *= scale;
    }
}

void svecd_print(svecd_t *a, const char *fmt)
{
    for (int i = 0; i < a->length; i++)
        printf(fmt, svecd_get(a, i));
    printf("\n");
}


 void svecd_add_inplace_keep(svecd_t *a, svecd_t *b, double bscale)
{

    int maxsz = a->nz + b->nz;
    // XXX beware of stack overflow
    TYPE tmp_values[maxsz];
    int tmp_indices[maxsz];
    int tmp_nz = 0;

    int aidx = 0;
    int bidx = 0;

    int ai = (aidx < a->nz) ? a->indices[aidx] : INT32_MAX;
    int bi = (bidx < b->nz) ? b->indices[bidx] : INT32_MAX;

    int start = 0;
    while (aidx < a->nz) {
        if(ai >= bi) {
            break;
        }
        aidx++;
        ai = (aidx < a->nz) ? a->indices[aidx] : INT32_MAX;
        start++;
    }

    while (aidx < a->nz || bidx < b->nz) {
        if (ai == bi) {
            // add and increment both.
            double val = a->values[aidx] + b->values[bidx] * bscale;
            if(fabs(val) > 1E-15) {
                tmp_values[tmp_nz] = val;
                tmp_indices[tmp_nz] = ai;
                tmp_nz++;
            }
            aidx++;
            ai = (aidx < a->nz) ? a->indices[aidx] : INT32_MAX;
            bidx++;
            bi = (bidx < b->nz) ? b->indices[bidx] : INT32_MAX;
        } else if (ai < bi) {
            if(fabs(a->values[aidx]) > 1E-15) {
                tmp_values[tmp_nz] = a->values[aidx];
                tmp_indices[tmp_nz] = ai;
                tmp_nz++;
            }
            aidx++;
            ai = (aidx < a->nz) ? a->indices[aidx] : INT32_MAX;
        } else {
            // bi < ai
            double val = b->values[bidx] * bscale;
            if(fabs(val) > 1E-15) {
                tmp_values[tmp_nz] = val;
                tmp_indices[tmp_nz] = bi;
                tmp_nz++;
            }
            bidx++;
            bi = (bidx < b->nz) ? b->indices[bidx] : INT32_MAX;
        }
    }
    a->nz = start + tmp_nz;
    svecd_ensure_capacity(a, a->nz);

    memcpy(&(a->indices[start]), tmp_indices, sizeof(int) * tmp_nz);
    memcpy(&(a->values[start]), tmp_values, sizeof(TYPE) * tmp_nz);
}

 void svecd_add_i0_x(svecd_t *x, svecd_t *a, svecd_t *b, double bscale, int aidx, int bidx)
{
    int maxsz = a->nz - aidx + b->nz - bidx;

    // XXX beware of stack overflow
    TYPE tmp_values[maxsz];
    int tmp_indices[maxsz];
    int tmp_nz = 0;

    int ai = (aidx < a->nz) ? a->indices[aidx] : INT32_MAX;
    int bi = (bidx < b->nz) ? b->indices[bidx] : INT32_MAX;

    while (aidx < a->nz || bidx < b->nz) {
        if (ai == bi) {
            // add and increment both.
            tmp_values[tmp_nz] = a->values[aidx] + b->values[bidx] * bscale;
            tmp_indices[tmp_nz] = ai;
            tmp_nz++;
            aidx++;
            ai = (aidx < a->nz) ? a->indices[aidx] : INT32_MAX;
            bidx++;
            bi = (bidx < b->nz) ? b->indices[bidx] : INT32_MAX;
        } else if (ai < bi) {
            tmp_values[tmp_nz] = a->values[aidx];
            tmp_indices[tmp_nz] = ai;
            tmp_nz++;
            aidx++;
            ai = (aidx < a->nz) ? a->indices[aidx] : INT32_MAX;
        } else {
            // bi < ai
            tmp_values[tmp_nz] = b->values[bidx] * bscale;
            tmp_indices[tmp_nz] = bi;
            tmp_nz++;
            bidx++;
            bi = (bidx < b->nz) ? b->indices[bidx] : INT32_MAX;
        }
    }

    svecd_ensure_capacity(x, tmp_nz);

    memcpy(x->indices, tmp_indices, sizeof(int)*tmp_nz);
    memcpy(x->values, tmp_values, sizeof(TYPE)*tmp_nz);
    x->nz = tmp_nz;
}

// Add vectors 'a' and 'b', treating elements in 'a' before index a0
// as zero, and elements in 'b' before index b0 as zero. Elements in b
// are scaled by 'bscale'.  In-place modifications are allowed (i.e.,
// x may be equal to a or b).
void svecd_add_i0(svecd_t *x, const svecd_t *a, const svecd_t *b, double bscale, int a0, int b0)
{
    int aidx = 0, bidx = 0;

    while (aidx < a->nz && a->indices[aidx] < a0)
        aidx++;
    while (bidx < b->nz && b->indices[bidx] < b0)
        bidx++;

    int maxsz = a->nz - aidx + b->nz - bidx;

    // XXX beware of stack overflow
    TYPE tmp_values[maxsz];
    int tmp_indices[maxsz];
    int tmp_nz = 0;

    // original code
    int ai, bi;

    int *aidxs = a->indices;
    int *bidxs = b->indices;

    double *avals = a->values;
    double *bvals = b->values;

    ai = (aidx < a->nz) ? aidxs[aidx] : INT32_MAX;
    bi = (bidx < b->nz) ? bidxs[bidx] : INT32_MAX;

    while (aidx < a->nz || bidx < b->nz) {
        if (ai == bi) {
            // add and increment both.
            tmp_values[tmp_nz] = avals[aidx] + bvals[bidx] * bscale;
            tmp_indices[tmp_nz] = ai;
            tmp_nz++;
            aidx++;
            ai = (aidx < a->nz) ? aidxs[aidx] : INT32_MAX;
            bidx++;
            bi = (bidx < b->nz) ? bidxs[bidx] : INT32_MAX;
        } else if (ai < bi) {
            tmp_values[tmp_nz] = avals[aidx];
            tmp_indices[tmp_nz] = ai;
            tmp_nz++;
            aidx++;
            ai = (aidx < a->nz) ? aidxs[aidx] : INT32_MAX;
        } else {
            // bi < ai
            tmp_values[tmp_nz] = bvals[bidx] * bscale;
            tmp_indices[tmp_nz] = bi;
            tmp_nz++;
            bidx++;
            bi = (bidx < b->nz) ? bidxs[bidx] : INT32_MAX;
        }
    }
    svecd_ensure_capacity(x, tmp_nz);

    memcpy(x->indices, tmp_indices, sizeof(int)*tmp_nz);
    memcpy(x->values, tmp_values, sizeof(TYPE)*tmp_nz);
    x->nz = tmp_nz;
}

/////////////////////////////////////////////////
// matrix operations

smatd_t *smatd_create(int nrows, int ncols)
{
    smatd_t *m = calloc(1, sizeof(smatd_t));
    m->nrows = nrows;
    m->ncols = ncols;

    // don't allocate anything for the rows yet.
    m->rows = calloc(nrows, sizeof(svecd_t));
    for (int i = 0; i < nrows; i++) {
        m->rows[i].length = ncols;
    }
    return m;
}

smatd_t *smatd_create_data(int nrows, int ncols, const TYPE *data)
{
    smatd_t *m = smatd_create(nrows, ncols);
    int pos = 0;

    for (int row = 0; row < nrows; row++) {
        for (int col = 0; col < ncols; col++) {
            smatd_set(m, row, col, data[pos]);
            pos++;
        }
    }

    return m;
}

smatd_t *smatd_copy(smatd_t *a)
{
    smatd_t *x = smatd_create(a->nrows, a->ncols);

    for (int i = 0; i < a->nrows; i++) {
        svecd_ensure_capacity(&x->rows[i], a->rows[i].nz);

        x->rows[i].nz = a->rows[i].nz;
        memcpy(x->rows[i].indices, a->rows[i].indices, sizeof(int)*a->rows[i].nz);
        memcpy(x->rows[i].values, a->rows[i].values, sizeof(TYPE)*a->rows[i].nz);
    }

    return x;
}

smatd_t *smatd_identity(int nrows, int ncols)
{
    smatd_t *m = smatd_create(nrows, ncols);

    int mx = nrows < ncols ? nrows : ncols;

    for (int i = 0; i < mx; i++)
        smatd_set(m, i, i, 1);

    return m;
}

void smatd_destroy(smatd_t *m)
{
    if (m == NULL)
        return;

    for (int i = 0; i < m->nrows; i++) {
        svecd_t *row = &m->rows[i];
        if (row) {
            free(row->indices);
            free(row->values);
            // don't free row; it's allocated as an array
        }
    }
    free(m->rows);
    free(m);
}

void smatd_set(smatd_t *m, int row, int col, TYPE val)
{
    svecd_set(&m->rows[row], col, val);
}

TYPE smatd_get(smatd_t *m, int row, int col)
{
    return svecd_get(&m->rows[row], col);
}

void smatd_print(smatd_t *m, const char *fmt)
{
    for (int i = 0; i < m->nrows; i++) {
        for (int j = 0; j < m->ncols; j++) {
            printf(fmt, smatd_get(m, i, j));
        }
        printf("\n");
    }
}

void smatd_debug(smatd_t *m, const char *fmt)
{
    printf("Matrix %d x %d\n", m->nrows, m->ncols);

    for (int i = 0; i < m->nrows; i++) {
        svecd_t *row = &m->rows[i];
        assert(row->length == m->ncols);

        for (int pos = 0; pos < row->nz; pos++) {
            printf("(%d: %d ", pos, row->indices[pos]);
            printf(fmt, row->values[pos]);
            printf(") ");
        }
        printf("\n");
    }
}

// CAUTION: Returns the internal representation of this
// matrix. Modifications of this vector will affect the matrix. Do not
// attempt to destroy the underlying vector.
svecd_t *smatd_get_row_volatile(smatd_t *m, int row)
{
    return &m->rows[row];
}

// unlike get_row, this returns a copy of the underlying data. This
// vector must be destroyed by the caller.
svecd_t *smatd_copy_column(smatd_t *m, int column)
{
    svecd_t *v = svecd_create(m->nrows);

    for (int row = 0; row < m->nrows; row++) {
        svecd_set(v, row, svecd_get(&m->rows[row], column));
    }

    return v;
}

// multiply using accumulation. Advantage is that we can skip
// completely over many zero cells. Disadvantage is that we skip
// around rows of b quite a lot and must accumulate the result (which
// means we have to work at keeping CSR order.)
//
// ---a--- ---x---  = a0*x + a1*y + a2*z
// ---b--- ---y---  = b0*x + b1*y + b2*z
// ---c--- ---z---  = c0*x + c1*y + c2*z
//
// On highly sparse matrices, this is much faster than the dotproduct
// method. On dense matrices, it's not a lot slower ==> make this the
// default.
smatd_t *smatd_multiply(smatd_t *a, smatd_t *b)
{
    assert(a->ncols == b->nrows);

    int xrows = a->nrows;
    int xcols = b->ncols;

    smatd_t *x = smatd_create(xrows, xcols);

    for (int i = 0; i < a->nrows; i++) {
        svecd_t *arow = &a->rows[i];
        svecd_t *xrow = &x->rows[i];

        // QQQ: Could also accumulate xrow using a dense vector, then collapse at the end.

        for (int pos = 0; pos < arow->nz; pos++) {
            TYPE v = arow->values[pos];
            int acol = arow->indices[pos];
            if (v==0)
                continue;

            assert(acol < b->nrows);
            svecd_t *brow = &b->rows[acol];
            svecd_add_i0(xrow, xrow, brow, v, 0, 0);
        }
    }

    return x;
}

// multiply using dot products.
smatd_t *smatd_multiply_dpmethod(smatd_t *a, smatd_t *b)
{
    assert(a->ncols == b->nrows);

    int xrows = a->nrows;
    int xcols = b->ncols;

    smatd_t *x = smatd_create(xrows, xcols);
    for (int col = 0; col < xcols; col++) {
        svecd_t *bcol = smatd_copy_column(b, col);

        for (int row = 0; row < xrows; row++) {
            svecd_t *arow = &a->rows[row];

            smatd_set(x, row, col, svecd_dot_product(arow, bcol));
        }

        svecd_destroy(bcol);
    }

    return x;
}

smatd_t *smatd_transpose(const smatd_t *a)
{
    smatd_t *x = smatd_create(a->ncols, a->nrows);
    int cursors[x->nrows];
    memset(cursors, 0, sizeof(cursors));

    for (int inrow = 0; inrow < a->nrows; inrow++) {
        svecd_t *arow = &a->rows[inrow];
        for (int pos = 0; pos < arow->nz; pos++) {
            int incol = arow->indices[pos];

            svecd_t *xrow = &x->rows[incol];
            svecd_ensure_capacity(xrow, ++cursors[incol]);
            xrow->indices[xrow->nz] = inrow;
            xrow->values[xrow->nz] = arow->values[pos];
            xrow->nz++;
        }
    }

    return x;
}

smatd_t *smatd_upper_right(smatd_t *a)
{
    smatd_t *x = smatd_create(a->nrows, a->ncols);
    for (int row = 0; row < a->nrows; row++) {
        svecd_t *av = &a->rows[row];
        svecd_t *xv = &x->rows[row];

        for (int aidx = 0; aidx < av->nz; aidx++) {
            if (av->indices[aidx] >= row) {
                xv->length = a->ncols;
                svecd_ensure_capacity(xv, av->nz - aidx);
                xv->nz = av->nz - aidx;
                memcpy(xv->indices, &av->indices[aidx], xv->nz * sizeof(int));
                memcpy(xv->values, &av->values[aidx], xv->nz * sizeof(TYPE));
                break;
            }
        }
    }
    return x;
}

// explanation of Cholesky approach:
//
// Consider matrix M (which we wish to factor).
//
// M =  [ A  B ]   Where A is 1x1, B is 1xN, C is NxN.
//      [ B' C ]
//
// The Cholesky factorization for this must be of the form:
//
// [ X  0 ] [ X  Y ] (where X is 1x1, Y is 1xN, etc.)
// [ Y' Z'] [ 0  Z ]
//
//  We must have X*X = A, or X = sqrt(A). Similarly, we must have Y =
//  B/sqrt(A). These choices suffice to completely (and correctly)
//  result in upper left portion (A, B, B') of M.
//
// But these choices for X and Y also contribute terms to the C parts
// of matrix M (due to the product of Y' and Y). So we subtract these
// contributions from C and then recursively factor what's left over.
//
// Specifically, C = Y'Y + Z'Z, so Z'Z = C - Y'Y. (We will then
// recursively factor Z'Z.)
smatd_chol_t *smatd_chol(smatd_t *a)
{
    // create a new matrix that is the upper-right of A.
    smatd_t *u = smatd_upper_right(a);

    int is_spd = (a->nrows == a->ncols);

    // recursively factor each "outer" upper-left shell (X,Y,Z),
    // modifying the internal bits.
    for (int i = 0; i < a->nrows; i++) {
        svecd_t *urowi = &u->rows[i];

        double d = svecd_get(urowi, i);
        is_spd &= (d > 0);
        d = sqrt(d);

        svecd_scale(urowi, 1.0 / d);

        // iterate over the non-zero elements... each non-zero element
        // will contribute something to block C, which we need to subtract.
        for (int pos = 0; pos < urowi->nz; pos++) {
            int j = urowi->indices[pos];
            if (j > i) {
                // apply an update to a row below us
                double s = urowi->values[pos];

                svecd_t *urowj = &u->rows[j];
//                assert(urowj->indices[0] == j);
//                assert(urowi->indices[pos] == j);
                svecd_add_i0_x(urowj, urowj, urowi, -s, 0, pos);
            }
        }
    }

    smatd_chol_t *chol = calloc(1, sizeof(smatd_chol_t));
    chol->is_spd = is_spd;
    chol->u = u;
    return chol;
}

void smatd_chol_reconstruct(smatd_chol_t *chol, int start_row)
{
    //reconstruct partial A?
    smatd_t *u = chol->u;
    for (int i = u->nrows - 1; i >= start_row; i--) {

        svecd_t *urowi = &u->rows[i];
        double d = svecd_get(urowi, i);

        // iterate over the non-zero elements... each non-zero element
        // will contribute something to block C, which we need to add it back.
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
}

void smatd_chol_inc(smatd_chol_t *chol, int start_row)
{
    smatd_t *u = smatd_upper_right(chol->u);
    chol->u = u;
    //smatd_t *u = chol->u;
    int is_spd = 1;

    // recursively factor each "outer" upper-left shell (X,Y,Z),
    // modifying the internal bits.
    for (int i = start_row; i < u->nrows; i++) {
        svecd_t *urowi = &u->rows[i];

        double d = svecd_get(urowi, i);
        is_spd &= (d > 0);
        d = sqrt(d);

        svecd_scale(urowi, 1.0 / d);
        // iterate over the non-zero elements... each non-zero element
        // will contribute something to block C, which we need to subtract.
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
    chol->is_spd = is_spd;
}

void smatd_chol_destroy(smatd_chol_t *chol)
{
    smatd_destroy(chol->u);
    free(chol);
}

smatd_qr_t *smatd_qr(const smatd_t *A, const double *_b)
{
    timeprofile_t *tp = timeprofile_create();
    timeprofile_stamp(tp, "begin");
    smatd_qr_t *qr = calloc(1, sizeof(smatd_qr_t));
    smatd_t *R = smatd_transpose(A);
    double *b = calloc(A->nrows, sizeof(double));
    memcpy(b, _b, A->nrows*sizeof(double));
    timeprofile_stamp(tp, "transpose");
//TODO: Given rotations

    //HouseHolder
    //R is m by n
    for(int k = 0; k < R->nrows; k++) {
        //select col k and row k:m
        svecd_t *v = svecd_create(R->ncols - k);
        svecd_t *rowk = &R->rows[k];

        //copy elements from rowk to v
        int nz_start_idx = 0;
        for (int i = 0; i < rowk->nz; i++) {
            if(rowk->indices[i] >= k) {
                nz_start_idx = i;
                break;
            }
        }

        if(rowk->indices[nz_start_idx] == k) {
            //(k,k) element is nonzero
            int v_nz = rowk->nz - nz_start_idx;
            svecd_ensure_capacity(v, v_nz);
            v->nz = v_nz;
            memcpy(v->indices, &rowk->indices[nz_start_idx], v_nz*sizeof(int));
            memcpy(v->values, &rowk->values[nz_start_idx], v_nz*sizeof(double));
        } else {
            //(k,k) element is zero
            int v_nz = rowk->nz - nz_start_idx;
            svecd_ensure_capacity(v, v_nz + 1);
            v->nz = v_nz + 1;
            v->indices[0] = k;
            memcpy(&v->indices[1], &rowk->indices[nz_start_idx], v_nz * sizeof(int));
            memcpy(&v->values[1], &rowk->values[nz_start_idx], v_nz * sizeof(double));
        }
        double sign = 1;
        double x_mag = 0;

        if(v->values[0] < 0)
            sign = -1;
        double v_mag_s = 0;
        for (int i = 0; i < v->nz; i++) {
            v_mag_s += v->values[i] * v->values[i];
        }
        double v_mag = sign * sqrt(v_mag_s);
        x_mag = v_mag;
        double v_0 = v->values[0] + v_mag;
        v_mag_s = v_mag_s - v->values[0]*v->values[0] + v_0*v_0;
        v_mag = sqrt(v_mag_s);
        v->values[0] = v_0;
        svecd_scale(v, 1.0 / v_mag);
        //when col_i = k; R(k,k) = -sign*v_mag; the others are all 0s;
        rowk->indices[nz_start_idx] = k;
        rowk->values[nz_start_idx] = -x_mag;
        rowk->nz = nz_start_idx + 1;
        //printf("rowk nz: %d, v nz:%d, k:%d, rows:%d, cols: %d\n", rowk->nz, v->nz, k, R->nrows, R->ncols);
        //printf("\n");
        for (int row_i = k+1; row_i < R->nrows; row_i++) {
            double vscale = 2*svecd_dot_product(v, &R->rows[row_i]);
            svecd_add_inplace_keep(&R->rows[row_i], v, -vscale);
        }
        double bscale = 0;
        for (int nz = 0; nz < v->nz; nz++) {
            int idx = v->indices[nz];
            bscale += 2*v->values[nz] * b[idx];
        }
        for (int nz = 0; nz < v->nz; nz++) {
            int idx = v->indices[nz];
            b[idx] = b[idx] - bscale * v->values[nz];
        }
        svecd_destroy(v);

    }
    timeprofile_stamp(tp, "find R");
    qr->R = smatd_transpose(R);
    timeprofile_stamp(tp, "transpose R");
    smatd_destroy(R);
    qr->b = b;
    //timeprofile_display(tp);
    timeprofile_destroy(tp);
    return qr;
}

void smatd_qr_solve(const smatd_qr_t *qr, double *x)
{
    //qr->R * x = qr->b;
    smatd_t *u = qr->R;
    double *b = qr->b;
    for (int i = u->ncols-1; i >= 0; i--) {
        double bi = b[i];
        svecd_t *urowi = &u->rows[i];
        assert(urowi->nz > 0);

        int jmp = 0;
        for (int idx = 0; idx<urowi->nz; idx++) {
            if(urowi->indices[idx] == i)
                break;
            else
                jmp++;
        }
        //TODO: need solve some numerical issues;
        //assert(urowi->indices[0] == i);

        TYPE diag = urowi->values[0+jmp];

        for (int pos = 1+jmp; pos < urowi->nz ; pos++)
            bi -= urowi->values[pos]*x[urowi->indices[pos]];

        x[i] = bi / diag;
    }
}

void smatd_qr_inc(smatd_qr_t *qr, smatd_t *meas, double *_b)
{
    smatd_t *R = qr->R;
    double *b = qr->b;
    assert(meas->ncols >= R->ncols);
    timeprofile_t *tp = timeprofile_create();
    timeprofile_stamp(tp, "begin");
    smatd_t *new_R = smatd_create(R->nrows + meas->nrows, meas->ncols);
    double *new_b = calloc(R->nrows + meas->nrows, sizeof(double));

    int start_row = INT32_MAX;

    for(int row_i = 0; row_i < R->nrows; row_i++) {
        if(R->rows[row_i].nz == 0)
            continue;
        svecd_t *row = &new_R->rows[row_i];
        svecd_ensure_capacity(row, R->rows[row_i].nz);
        row->nz = R->rows[row_i].nz;
        memcpy(row->indices, R->rows[row_i].indices, sizeof(int)*R->rows[row_i].nz);
        memcpy(row->values, R->rows[row_i].values, sizeof(double)*R->rows[row_i].nz);
    }
    for(int row_i = 0; row_i < meas->nrows; row_i++) {
        if(meas->rows[row_i].nz == 0)
            continue;
        if(start_row > meas->rows[row_i].indices[0])
            start_row = meas->rows[row_i].indices[0];

        svecd_t *row = &new_R->rows[row_i + R->nrows];
        svecd_ensure_capacity(row, meas->rows[row_i].nz);
        row->nz = meas->rows[row_i].nz;
        memcpy(row->indices, meas->rows[row_i].indices, sizeof(int) * meas->rows[row_i].nz);
        memcpy(row->values, meas->rows[row_i].values, sizeof(double) * meas->rows[row_i].nz);
    }

    //new_b
    memcpy(new_b, b, sizeof(double) * R->nrows);
    memcpy(&new_b[R->nrows], _b, sizeof(double) * meas->nrows);

    smatd_destroy(R);
    free(b);
    R = smatd_transpose(new_R);
    smatd_destroy(new_R);
    b = new_b;

    timeprofile_stamp(tp, "new A");
    //printf("##########start-row:%d \n", start_row);

    for(int k = start_row; k < R->nrows; k++) {
        //select col k and row k:m
        svecd_t *v = svecd_create(R->ncols - k);
        svecd_t *rowk = &R->rows[k];

        //copy elements from rowk to v
        int nz_start_idx = 0;
        for (int i = 0; i < rowk->nz; i++) {
            if(rowk->indices[i] >= k) {
                nz_start_idx = i;
                break;
            }
        }
        if(rowk->indices[nz_start_idx] == k) {
            //(k,k) element is nonzero
            int v_nz = rowk->nz - nz_start_idx;
            svecd_ensure_capacity(v, v_nz);
            v->nz = v_nz;
            memcpy(v->indices, &rowk->indices[nz_start_idx], v_nz*sizeof(int));
            memcpy(v->values, &rowk->values[nz_start_idx], v_nz*sizeof(double));
        } else {
            //(k,k) element is zero
            int v_nz = rowk->nz - nz_start_idx;
            svecd_ensure_capacity(v, v_nz + 1);
            v->nz = v_nz + 1;
            v->indices[0] = k;
            memcpy(&v->indices[1], &rowk->indices[nz_start_idx], v_nz * sizeof(int));
            memcpy(&v->values[1], &rowk->values[nz_start_idx], v_nz * sizeof(double));
        }
        double sign = 1;
        double x_mag = 0;

        if(v->values[0] < 0)
            sign = -1;
        double v_mag_s = 0;
        for (int i = 0; i < v->nz; i++) {
            v_mag_s += v->values[i] * v->values[i];
        }
        double v_mag = sign * sqrt(v_mag_s);
        x_mag = v_mag;
        double v_0 = v->values[0] + v_mag;
        v_mag_s = v_mag_s - v->values[0]*v->values[0] + v_0*v_0;
        v_mag = sqrt(v_mag_s);
        v->values[0] = v_0;
        svecd_scale(v, 1.0 / v_mag);

        //when col_i = k; R(k,k) = -sign*v_mag; the others are all 0s;
        rowk->indices[nz_start_idx] = k;
        rowk->values[nz_start_idx] = -x_mag;
        rowk->nz = nz_start_idx + 1;

        //printf("v nz:%d, k:%d, rows:%d, cols: %d\n", v->nz, k, R->nrows, R->ncols);
        //printf("\n");
        (void)x_mag;
        for (int row_i = k+1; row_i < R->nrows; row_i++) {
            double vscale = 2*svecd_dot_product(v, &R->rows[row_i]);
            svecd_add_inplace_keep(&R->rows[row_i], v, -vscale);
        }


        double bscale = 0;
        for (int nz = 0; nz < v->nz; nz++) {
            int idx = v->indices[nz];
            bscale += 2*v->values[nz] * b[idx];
        }
        for (int nz = 0; nz < v->nz; nz++) {
            int idx = v->indices[nz];
            b[idx] = b[idx] - bscale * v->values[nz];
        }
        svecd_destroy(v);
    }
    timeprofile_stamp(tp, "solve R inc");
    qr->R = smatd_transpose(R);
    timeprofile_stamp(tp, "transpose");
    //timeprofile_display(tp);
    timeprofile_destroy(tp);
    smatd_destroy(R);
    qr->b = b;
}

void smatd_qr_destroy(smatd_qr_t *qr)
{
    if(!qr)
        return;

    if(qr->b)
        free(qr->b);
    smatd_destroy(qr->R);
}

void smatd_ltriangle_solve(smatd_t *L, const TYPE *b, TYPE *x)
{
    for (int i = 0; i < L->ncols; i++) {
        double bi = b[i];

        svecd_t *Lrowi = &L->rows[i];

        // this row should end with a value along the diagonal.
        // NB: rank-deficient matrices could have a zero here?
        assert(Lrowi->nz > 0);
        assert(Lrowi->indices[Lrowi->nz-1] == i);

        TYPE diag = Lrowi->values[Lrowi->nz-1];

        // compute dot product of Lrowi and the first (i-1) elements
        // of y.  Note: exploit fact that l is exactly
        // triangular. Don't need range checks.
        for (int pos = 0; pos < Lrowi->nz - 1; pos++)
            bi -= Lrowi->values[pos]*x[Lrowi->indices[pos]];

        x[i] = bi / diag;
    }
}

// Solve the triangular system Lx = b, where L is provided in the form
// of its transpose (U). The math for this follows from solving the
// system: x'U = b' and proceeds recursively (similar to the cholesky
// factorization).
void smatd_ltransposetriangle_solve(smatd_t *u, const TYPE *b, TYPE *x)
{
    int n = u->ncols;

    memcpy(x, b, n*sizeof(TYPE));

    for (int i = 0; i < n; i++) {
        svecd_t *urowi = &u->rows[i];
        /* if(i == 0) { */
        /*     for (int i = 0; i < n; i++) { */
        /*         printf("%f,", x[i]); */
        /*     } */
        /*     puts("\n x[i]-----"); */
        /* } */
        x[i] = x[i] / urowi->values[0];

        for (int pos = 1; pos < urowi->nz; pos++) {
            int uidx = urowi->indices[pos];
            x[uidx] -= x[i] * urowi->values[pos];
        }
    }

}

void smatd_utriangle_solve(smatd_t *u, const TYPE *b, TYPE *x)
{
    for (int i = u->ncols-1; i >= 0; i--) {
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
    }
}

// Solve the system Ax=b. User provides storage for x.
void smatd_chol_solve_full(smatd_chol_t *chol, const TYPE *b, TYPE *y, TYPE *x)
{
    smatd_t *u = chol->u;

    // back solve U'y = b
    // smatd_t *l = smatd_transpose(u);
    // smatd_ltriangle_solve(l, b, y);
    // smatd_destroy(l);

    // this replacement is around 10x faster than transposing and using ltriangle_solve
    smatd_ltransposetriangle_solve(u, b, y);

    // back solve Ux = y
    smatd_utriangle_solve(u, y, x);
}

void smatd_chol_solve(smatd_chol_t *chol, const TYPE *b, TYPE *x)
{
    // U'Ux = b
    // U'y = b

    smatd_t *u = chol->u;

    TYPE y[u->ncols];

        // back solve U'y = b
        // smatd_t *l = smatd_transpose(u);
        // smatd_ltriangle_solve(l, b, y);
        // smatd_destroy(l);

    // this replacement is around 10x faster than transposing and using ltriangle_solve
    smatd_ltransposetriangle_solve(u, b, y);

    // back solve Ux = y
    smatd_utriangle_solve(u, y, x);
}

// basically the same idea as chol, except we don't assume symmetry.
//
// Begin with a block decomposition of M, with A = 1x1, B 1xN, C Nx1, D NxN.
// Consider just the LDU decomposition of this block decomposition:
//
// M = [ A B ] = [ I  0  ] [ X  0 ] [ I  Z ] = [ X       XZ   ]
//     [ C D ]   [ W  I  ] [ 0  Y ] [ 0  I ]   [ WX   WXZ + Y ]
//
// X = A
// Z = inv(A)*B
// W = C*inv(A)
// Y = D - C*inv(A)*B
//
// Or, in other words, we want to transform like this:
//
// [ A B ] => [ X Z ] = [ A            inv(A)*B     ]
// [ C D ]    [ W Y ]   [ C*inv(A)   D - C*inv(A)*B ]
//
// We then recurse on the Y submatrix.
//
// The result will be a single matrix with the L, D, and U all packed
// together (I and 0 terms omitted.)
smatd_ldu_t *smatd_ldu(smatd_t *a)
{
    smatd_t *ldu = smatd_copy(a);

    // recursively do LDU decomposition for submatrix of M with upper
    // corner at i,i:
    for (int i = 0; i < ldu->nrows; i++) {
        svecd_t *ldurowi = &ldu->rows[i];

        // get A block
        TYPE d = svecd_get(ldurowi, i);
        TYPE invd = 1.0 / d;

        // handle Z sub-block: the rest of the row is scaled by 1/d
        svecd_scale_i0(ldurowi, invd, i+1);

        for (int j = i+1; j < ldu->nrows; j++) {
            // handle W sub-block
            // INEFFICIENT?
            TYPE v = svecd_get(&ldu->rows[j], i);
            svecd_set(&ldu->rows[j], i, v * invd);

            // handle Y sub-block--- subtract C*inv(A)*B
            if (v == 0)
                continue;

            // Note that we've already scaled this ldurowi by invd.
            svecd_add_i0(&ldu->rows[j], &ldu->rows[j], ldurowi, -v, 0, i+1);
        }
    }

    smatd_ldu_t *sldu = calloc(1, sizeof(smatd_ldu_t));
    sldu->ldu = ldu;

    return sldu;
}

void smatd_ldu_destroy(smatd_ldu_t *sldu)
{
    if (sldu==NULL)
        return;
    smatd_destroy(sldu->ldu);
    free(sldu);
}

void smatd_ldu_get(smatd_ldu_t *sldu, smatd_t **l_out, smatd_t **d_out, smatd_t **u_out)
{
    smatd_t *ldu = sldu->ldu;

    int R = ldu->nrows < ldu->ncols ? ldu->nrows : ldu->ncols;

    smatd_t *l = smatd_create(ldu->nrows, R);
    smatd_t *d = smatd_create(R, R);
    smatd_t *u = smatd_create(R, ldu->ncols);

    for (int row = 0; row < ldu->nrows; row++) {
        svecd_t *ldurow = &ldu->rows[row];

        int lnz = 0;
        int unz = 0;
        int upos = 0;

        for (int pos = 0; pos < ldurow->nz; pos++) {

            if (ldurow->indices[pos] == row) {
                // diagonal element? add to D.
                if (row < R) {
                    svecd_t *drow = &d->rows[row];

                    svecd_ensure_capacity(drow, 1);
                    drow->nz = 1;
                    drow->indices[0] = row;
                    drow->values[0] = ldurow->values[pos];
                }
                continue;
            } else if (ldurow->indices[pos] < row) {
                lnz++;
                continue;
            } else {
                // index must be into the U matrix (indices[pos] > row)
                upos = pos;
                unz = ldurow->nz - upos;
                break; // don't need to do any more searching.
            }
        }

        if (row < l->nrows) {
            svecd_t *lrow = &l->rows[row];
            svecd_ensure_capacity(lrow, lnz+1);
            memcpy(&lrow->indices[0], &ldurow->indices[0], lnz*sizeof(int));
            memcpy(&lrow->values[0], &ldurow->values[0], lnz*sizeof(TYPE));
            if (row < lrow->length) {
                lrow->indices[lnz] = row;
                lrow->values[lnz] = 1;
                lrow->nz = lnz+1;
            }
        }

        if (row < u->nrows) {
            svecd_t *urow = &u->rows[row];
            svecd_ensure_capacity(urow, unz+1);
            memcpy(&urow->indices[1], &ldurow->indices[upos], unz*sizeof(int));
            memcpy(&urow->values[1], &ldurow->values[upos], unz*sizeof(TYPE));
            if (row < urow->length) {
                urow->indices[0] = row;
                urow->values[0] = 1;
                urow->nz = unz+1;
            }
        }
    }

    *l_out = l;
    *d_out = d;
    *u_out = u;
}

int smatd_nz(smatd_t *A)
{
    int nz = 0;

    for (int rowi = 0; rowi < A->nrows; rowi++) {
        svecd_t *row = &A->rows[rowi];

        nz += row->nz;
    }

    return nz;
}
