#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "AUX/LinearAlgebra.hpp"

namespace fatrop
{
    // cpy elements form sx to sy but in reversed order to avoid aliasing issues in recursion
    void fatrop_dveccp_reversed(int m, struct blasfeo_dvec *sx, int xi, struct blasfeo_dvec *sy, int yi)
    {
        for (int i = m - 1; i >= 0; i--)
        {
            VECEL(sy, yi + i) = VECEL(sx, xi + i);
        }
    }

    // void fatrop_potrf_l_mn(int m, int n, struct blasfeo_dmat *sC, int ci, int cj, struct blasfeo_dmat *sD, int di, int dj)
    // {
    //     blasfeo_dpotrf_l_mn(m, n, sC, ci, cj, sD, di, dj);
    //     int minmn = (m < n) ? m : n;
    //     for (int i =0; i<minmn; i++){
    //         assert(MATEL(sD, di+i,dj+i)>0);
    //     }
    // }
    /** \brief D <= alpha * B * A^{-1} , with A lower triangular employing explicit inverse of diagonal, fatrop uses its own (naive) implementation since it  not implemented yet in blasfeo*/
    void fatrop_dtrsm_rlnn(int m, int n, double alpha, MAT *sA, int offs_ai, int offs_aj, MAT *sB, int offs_bi, int offs_bj, MAT *sD, int offs_di, int offs_dj)
    {
        sD->use_dA = 0;
        for (int aj = n - 1; aj >= 0; aj--)
        {
            double ajj = MATEL(sA, offs_ai + aj, aj + offs_aj);
            double inv_ajj = 1.0 / ajj;
            double scjj = alpha * inv_ajj;
            for (int k = 0; k < m; k++)
            {
                MATEL(sD, offs_di + k, offs_dj + aj) = scjj * MATEL(sB, offs_bi + k, offs_bj + aj);
            }
            for (int ai = aj + 1; ai < n; ai++)
            {
                double sc = -inv_ajj * MATEL(sA, offs_ai + ai, offs_aj + aj);
                for (int k = 0; k < m; k++)
                {
                    // this algorithm is "store bounded", the loops can be switched like in the alt version to make this more efficient
                    MATEL(sD, offs_di + k, offs_dj + aj) += sc * MATEL(sD, offs_di + k, offs_dj + ai);
                }
            }
        }
    }
    // /** \brief D <= alpha * B * A^{-1} , with A lower triangular employing explicit inverse of diagonal, fatrop uses its own (naive) implementation since it  not implemented yet in blasfeo*/
    // // this is an experimental, more efficient, but the corner cases are not treated in unrolled loop!! it achieves a speed-up of about factor 3 w.r.t naive implementation
    // void fatrop_dtrsm_rlnn_alt(int m, int n, double alpha, MAT *sA, int offs_ai, int offs_aj, MAT *sB, int offs_bi, int offs_bj, MAT *sD, int offs_di, int offs_dj)
    // {
    //     for (int aj = n - 1; aj >= 0; aj--)
    //     {
    //         double ajj = MATEL(sA, offs_ai + aj, aj + offs_aj);
    //         double inv_ajj = 1.0 / ajj;
    //         double scjj = alpha * inv_ajj;
    //         for (int k = 0; k < m; k++)
    //         {
    //             // todo, check if possible to incude in main loop
    //             MATEL(sD, offs_di + k, offs_dj + aj) = scjj * MATEL(sB, offs_bi + k, offs_bj + aj);
    //         }
    //         for (int k = 0; k < m; k++)
    //         {
    //             double res = 0.0;
    //             double res1 = 0.0;
    //             double res2 = 0.0;
    //             double res3 = 0.0;
    //             // double res4 = 0.0;
    //             // double res5 = 0.0;
    //             // double res6 = 0.0;
    //             // double res7 = 0.0;
    //             for (int ai = aj + 1; ai < n; ai = ai + 4)
    //             {
    //                 // todo unroll loop -> more independent operations -> filled pipelines
    //                 double sc = -inv_ajj * MATEL(sA, offs_ai + ai, offs_aj + aj);
    //                 res += sc * MATEL(sD, offs_di + k, offs_dj + ai);
    //                 double sc1 = -inv_ajj * MATEL(sA, offs_ai + ai + 1, offs_aj + aj);
    //                 res1 += sc1 * MATEL(sD, offs_di + k, offs_dj + ai + 1);
    //                 double sc2 = -inv_ajj * MATEL(sA, offs_ai + ai + 2, offs_aj + aj);
    //                 res2 += sc2 * MATEL(sD, offs_di + k, offs_dj + ai + 2);
    //                 double sc3 = -inv_ajj * MATEL(sA, offs_ai + ai + 3, offs_aj + aj);
    //                 res3 += sc3 * MATEL(sD, offs_di + k, offs_dj + ai + 3);
    //                 // double sc4 = -inv_ajj * MATEL(sA, offs_ai + ai+4, offs_aj + aj);
    //                 // res4 += sc * MATEL(sD, offs_di + k, offs_dj + ai+4);
    //                 // double sc5 = -inv_ajj * MATEL(sA, offs_ai + ai+5, offs_aj + aj);
    //                 // res5 += sc * MATEL(sD, offs_di + k, offs_dj + ai+5);
    //                 // double sc6 = -inv_ajj * MATEL(sA, offs_ai + ai+6, offs_aj + aj);
    //                 // res6 += sc * MATEL(sD, offs_di + k, offs_dj + ai+6);
    //                 // double sc7 = -inv_ajj * MATEL(sA, offs_ai + ai+7, offs_aj + aj);
    //                 // res7 += sc * MATEL(sD, offs_di + k, offs_dj + ai+7);
    //             }
    //             // MATEL(sD, offs_di + k, offs_dj + aj) += (res+res1+res2+res3+res4+res5+res6+res7);
    //             MATEL(sD, offs_di + k, offs_dj + aj) += (res + res1) + (res2 + res3);
    //         }
    //     }
    // }
    // B <= B + alpha*A^T (B is mxn)
    void fatrop_dgead_transposed(int m, int n, double alpha, struct blasfeo_dmat *sA, int offs_ai, int offs_aj, struct blasfeo_dmat *sB, int offs_bi, int offs_bj)
    {
        for (int bj = 0; bj < n; bj++)
        {
            for (int bi = 0; bi < m; bi++)
            {
                MATEL(sB, offs_bi + bi, offs_bj + bj) += alpha * MATEL(sA, offs_ai + bj, offs_aj + bi);
            }
        }
    }

    /** \brief returns the maximum element of a blasfeo matrix of size (m,n), starting at (ai,aj) */
    MatrixInd max_el(int m, int n, MAT *matr, int ai, int aj)
    {
        MatrixInd res;
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
    /** \brief Function to calculate LU factorization result is saved in A, L is lower unitriangular */
    void LU_FACT(const int m, const int n, const int n_max, int &rank, MAT *A, PMAT *Pl_p, PMAT *Pr_p, double tol)
    {
        A->use_dA = 0;
        int *perm_left = (int *)(*Pl_p);
        int *perm_right = (int *)(*Pr_p);
        int minmn = MIN(m, n_max);
        int j = 0;
        for (int i = 0; i < minmn; i++)
        {
            MatrixInd max_curr = max_el(m, n_max, A, i, i);
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
        At->use_dA = 0;
        int *perm_left = (int *)(*Pl_p);
        int *perm_right = (int *)(*Pr_p);
        int minmn = MIN(m, n_max);
        int j = 0;
        for (int i = 0; i < minmn; i++)
        {
            MatrixInd max_curr = max_el(n_max, m, At, i, i);
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

    void fatrop_dtrsv_unu(const int m, const int n, blasfeo_dmat *sA, const int ai, const int aj, blasfeo_dvec *sx, const int xi, blasfeo_dvec *sz, const int zi)
    {
        for (int i = m; i < n; i++)
        {
            VECEL(sz, zi + i) = VECEL(sx, xi + i);
        }
        for (int i = m - 1; i >= 0; i--)
        {
            double res = VECEL(sx, xi + i);
            for (int j = i + 1; j < n; j++)
            {
                res -= MATEL(sA, ai + i, aj + j) * VECEL(sz, zi + j);
            }
            VECEL(sz, zi + i) = res;
        }
    }
    void fatrop_identity(const int m, MAT *sA, const int ai, const int aj)
    {
        GESE(m, m, 0.0, sA, ai, aj);
        DIARE(m, 1.0, sA, ai, aj);
    }
    void fatrop_drowad(int kmax, double alpha, struct blasfeo_dvec *sx, int xi, struct blasfeo_dmat *sA, int ai, int aj)
    {
        for (int i = 0; i < kmax; i++)
        {
            MATEL(sA, ai, aj + i) += alpha* VECEL(sx, xi + i);
        }
    }
} // namespace fatrop