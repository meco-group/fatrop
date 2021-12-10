#include "FatropLinearAlgebraBlasfeo.hpp"


namespace fatrop{

    /** \brief D <= alpha * B * A^{-1} , with A lower triangular employing explicit inverse of diagonal, fatrop uses its own (naive) implementation since it  not implemented yet in blasfeo*/
    void fatrop_dtrsm_rlnn(int m, int n, double alpha, MAT *sA, int offs_ai, int offs_aj, MAT *sB, int offs_bi, int offs_bj, MAT *sD, int offs_di, int offs_dj)
    {
        for (int aj = n - 1; aj >= 0; aj--)
        {
            double ajj = MATEL(sA, aj + offs_ai, aj + offs_aj);
            double sc = alpha / ajj;
            for (int k = 0; k < m; k++)
            {
                MATEL(sD, k + offs_di, aj + offs_dj) = sc * MATEL(sB, k + offs_bi, aj + offs_bj);
            }
            for (int ai = aj + 1; ai < m; ai++)
            {
                double sc = -alpha * MATEL(sA, ai + offs_ai, aj + offs_aj) / ajj;
                for (int k = 0; k < m; k++)
                {
                    MATEL(sD, k + offs_di, aj + offs_dj) += sc * MATEL(sD, k + offs_di, ai + offs_dj);
                }
            }
        }
    }
    /** \brief returns the maximum element of a blasfeo matrix of size (m,n), starting at (ai,aj) */
    matrix_ind max_el(int m, int n, MAT *matr, int ai, int aj)
    {
        matrix_ind res;
        res.ai = ai;
        res.aj = aj;
        double valmax = 0.0;
        for (int j = aj; j < n; j++)
        {
            for (int i = ai; i < m; i++)
            {
                double valij = abs(MATEL(matr, i, j));
                if (valij >= valmax)
                {
                    valmax = valij;
                    res.ai = i;
                    res.aj = j;
                }
            }
        }
        return res;
    };
    /** \brief Function to calculate LU factorization result is saved in A, L is unit diagonal */
    void LU_FACT(const int m, const int n, const int n_max, int &rank, MAT *A, PMAT *Pl_p, PMAT *Pr_p, double tol)
    {
        int *perm_left = (int *)(*Pl_p);
        int *perm_right = (int *)(*Pr_p);
        int minmn = MIN(m, n_max);
        int j = 0;
        for (int i = 0; i < minmn; i++)
        {
            matrix_ind max_curr = max_el(m, n_max, A, i, i);
            if (abs(MATEL(A, max_curr.ai, max_curr.aj)) < tol)
            {
                break;
            }
            // switch rows
            ROWSW(n, A, i, 0, A, max_curr.ai, 0);
            // save in permutation vector
            perm_left[i] = max_curr.ai;
            // switch cols
            COLSW(m, A, 0, i, A, 0, max_curr.aj);
            // save in permutation vector
            perm_right[i] = max_curr.aj;
            for (int j = i + 1; j < m; j++)
            {
                double Lji = MATEL(A, j, i) / MATEL(A, i, i);
                MATEL(A, j, i) = Lji;
                GEAD(1, n - (i + 1), -Lji, A, i, i + 1, A, j, i + 1);
            }
            j = i + 1;
        }
        rank = j;
    };
    /** \brief Function to calculate LU factorization but A, and result (L and U) are transposed, all indices refer to the dimensions of the original A matrix (and not the transposed one) */
    void LU_FACT_transposed(const int m, const int n, const int n_max, int &rank, MAT *At, PMAT *Pl_p, PMAT *Pr_p, double tol)
    {
        int *perm_left = (int *)(*Pl_p);
        int *perm_right = (int *)(*Pr_p);
        int minmn = MIN(m, n_max);
        int j = 0;
        for (int i = 0; i < minmn; i++)
        {
            matrix_ind max_curr = max_el(n_max, m, At, i, i);
            if (abs(MATEL(At, max_curr.ai, max_curr.aj)) < tol)
            {
                break;
            }
            // switch rows
            COLSW(n, At, 0, i, At, 0, max_curr.aj);
            // save in permutation vector
            perm_left[i] = max_curr.aj;
            // switch cols
            ROWSW(m, At, i, 0, At, max_curr.ai, 0);
            // save in permutation vector
            perm_right[i] = max_curr.ai;
            for (int j = i + 1; j < m; j++)
            {
                double Lji = MATEL(At, i, j) / MATEL(At, i, i);
                MATEL(At, i, j) = Lji;
                GEAD(n - (i + 1), 1, -Lji, At, i + 1, i, At, i + 1, j);
            }
            j = i + 1;
        }
        rank = j;
    };
} // namespace fatrop