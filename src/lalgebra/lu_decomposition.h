#ifndef LU_DECOMPOSITION_H
#define LU_DECOMPOSITION_H

#include "../dr.h"

    struct LUdecomposition {
        int n;
        DR::matrix_double lu;
        DR::vector_int index;
        double d;
        LUdecomposition(DR::matrix_double_I &a);
        void solve(DR::vector_double_I &b, DR::vector_double_O &x);
        void solve(DR::matrix_double_O & b, DR::matrix_double_O &x);
        void inverse(DR::matrix_double_O &ainv);
        double det();
        void mprove(DR::vector_double_O &b, DR::vector_double_IO &x);
        DR::matrix_double_I &aref;
    };

LUdecomposition::LUdecomposition(DR::matrix_double_I &a) :
    n(a.nrows()), lu(a), aref(a), index(n)
{
    const double TINY = 1.0e-40;
    int imax, i, j, k;
    double big, temp;
    DR::vector_double vv(n);
    d = 1.0;

    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++) {
            temp = abs(lu[i][j]);
            if (temp > big) big = temp;
        }
        if (big == 0.0) throw("Singular matrix in LUdecomposition");

        // Save the scaling
        vv[i] = 1.0 / big;
    }

    for (k = 0; k < n; k++) {
        big = 0.0;
        for (i = k; i < n; i++) {
            temp = vv[i] * abs(lu[i][k]);
            if (temp > big) {
                big = temp;
                imax = i;
            }
        }

        // Do we need to interchange rows? Yes, do so...
        if (k != imax) {
            for (j = 0; j < n; j++) {
                temp = lu[imax][j];
                lu[imax][j] = lu[k][j];
                lu[k][j] = temp;
            }
            d = -d;

            // Also interchange the scale factor
            vv[imax] = vv[k];
        }
        index[k] = imax;
        if (lu[k][k] == 0.0) lu[k][k] = TINY;

        for (i = k + 1; i < n; i++) {
            lu[i][k] /= lu[k][k];
            temp = lu[i][k];

            for (j = k + 1; j < n; j++)
                lu[i][j] -= temp * lu[k][j];
        }
    }
}

#endif // LU_DECOMPOSITION_H
