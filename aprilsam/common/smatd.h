/* $COPYRIGHT_UM
$LICENSE_LGPL
*/

/*$LICENSE*/

#ifndef _SMATD_H
#define _SMATD_H

// sparse matrix library.  ebolson 9/2013, adapted from
// april.jmat. Unlike jmat, doesn't support per-row variation of
// sparse/dense representation. However, API supports some specific
// dense operations with matd objects.

typedef struct
{
    int length; // logical length
    int alloc;  // allocated size of indices and values.

    int nz;     // non-zero length (and valid length of indices/values arrays)
    int *indices;
    double *values;

    int getcursor, setcursor;
} svecd_t;

typedef struct
{
    int nrows;
    int ncols;

    svecd_t *rows;
    svecd_t *cols;
} smatd_t;

#define TYPE double

smatd_t *smatd_create(int nrows, int ncols);
smatd_t *smatd_create_data(int nrows, int ncols, const TYPE *data);
smatd_t *smatd_identity(int nrows, int ncols);
void smatd_destroy(smatd_t *m);
void smatd_set(smatd_t *m, int row, int col, TYPE val);
TYPE smatd_get(smatd_t *m, int row, int col);

void smatd_print(smatd_t *m, const char *fmt);

int smatd_nz(smatd_t *A);

// CAUTION: Returns the internal representation of this
// matrix. Modifications of this vector will affect the matrix. Do not
// attempt to destroy the underlying vector.
svecd_t *smatd_get_row_volatile(smatd_t *m, int row);

// unlike get_row, this returns a copy of the underlying data. This
// vector must be destroyed by the caller.
svecd_t *smatd_copy_column(smatd_t *m, int column);

smatd_t *smatd_multiply(smatd_t *a, smatd_t *b);
smatd_t *smatd_upper_right(smatd_t *a);
smatd_t *smatd_transpose(const smatd_t *a);
void smatd_ltransposetriangle_solve(smatd_t *u, const TYPE *b, TYPE *x);
typedef struct
{
    smatd_t *u;
    int is_spd;

} smatd_chol_t;

typedef struct
{
    double *b;
    smatd_t *R;
} smatd_qr_t;

typedef struct
{
    smatd_t *ldu;
} smatd_ldu_t;

// Compute a Cholesky factorization. No pivoting is performed
// (deliberately--- so that we do not undermine an externally-applied
// sparsity-maximizing variable ordering. The matrix must be square and
// SPD, or the results are undefined.
smatd_chol_t *smatd_chol(smatd_t *a);
void smatd_utriangle_solve(smatd_t *u, const TYPE *b, TYPE *x);
void smatd_chol_solve(smatd_chol_t *chol, const TYPE *b, TYPE *x);
void smatd_chol_solve_full(smatd_chol_t *chol, const TYPE *b, TYPE *y, TYPE *x);
void smatd_chol_reconstruct(smatd_chol_t *chol, int start_row);
void smatd_chol_inc(smatd_chol_t *chol, int start_row);
void smatd_chol_destroy(smatd_chol_t *chol);

smatd_qr_t *smatd_qr(const smatd_t *A, const double *b);
void smatd_qr_inc(smatd_qr_t *qr, smatd_t *meas, double *_b);
void smatd_qr_solve(const smatd_qr_t *qr, double *x);
void smatd_qr_destroy(smatd_qr_t *qr);

// Compute an LDU factorization. No pivoting is performed
// (deliberately--- so that we do not undermine an externally-applied
// sparsity-maximizing variable ordering. The matrix can be any shape,
// and should have rank equal to its minimum dimension (e.g., a 10x3
// matrix should have rank 3.)
smatd_ldu_t *smatd_ldu(smatd_t *a);
// get the L, D, and U matrices computed by smatd_ldu. Note: not
// needed to solve systems.
// LDU exists for any matrix. L: MxR, D: RxR, U: RxN
// Where R is min(M,R)
void smatd_ldu_get(smatd_ldu_t *sldu, smatd_t **l_out, smatd_t **d_out, smatd_t **u_out);

void smatd_ldu_destroy(smatd_ldu_t *sldu);

void svecd_add_i0_x(svecd_t *x, svecd_t *a, svecd_t *b, double bscale, int aidx, int bidx);
TYPE svecd_get(svecd_t *v, int idx);
void svecd_scale(svecd_t *v, double scale);

#endif
