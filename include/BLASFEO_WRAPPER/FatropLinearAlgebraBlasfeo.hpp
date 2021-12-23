#ifndef FATROP_BLASFEO_INCLUDED
#define FATROP_BLASFEO_INCLUDED

// macros
#define MAT blasfeo_dmat
#define VEC blasfeo_dvec
#define MEMSIZE_MAT blasfeo_memsize_dmat
#define CREATE_MAT blasfeo_create_dmat
#define ROWPE blasfeo_drowpe
#define VECPE blasfeo_dvecpe
#define VECPEI blasfeo_dvecpei
#define ROWPEI blasfeo_drowpei
#define COLPE blasfeo_dcolpe
#define COLPEI blasfeo_dcolpei
#define CREATE_MAT blasfeo_create_dmat
#define MATEL BLASFEO_DMATEL
#define ROWSW blasfeo_drowsw
#define COLSW blasfeo_dcolsw
#define GEAD blasfeo_dgead
#define GEADTR fatrop_dgead_transposed
#define GECP blasfeo_dgecp
#define VECCP blasfeo_dveccp
#define TRSM_RLNN fatrop_dtrsm_rlnn //TODO this is not implemented by blasfeo so we defined our own (naive) implementation
#define VECEL BLASFEO_DVECEL
#define MEMSIZE_VEC blasfeo_memsize_dvec
#define CREATE_VEC blasfeo_create_dvec
#define GEMM_NT blasfeo_dgemm_nt
#define GEAD blasfeo_dgead
#define SYRK_LN_MN blasfeo_dsyrk_ln_mn
#define GETR blasfeo_dgetr
#define TRTR_L blasfeo_dtrtr_l
// #define POTRF_L_MN blasfeo_dpotrf_l_mn

#define POTRF_L_MN fatrop_potrf_l_mn
#define ROWEX blasfeo_drowex
#define TRSV_LTN blasfeo_dtrsv_ltn
#define TRSV_LNN blasfeo_dtrsv_lnn
#define TRSV_UNU fatrop_dtrsv_unu
#define GEMV_T blasfeo_dgemv_t
#define PMAT fatrop_permutation_matrix

#include <iostream>
extern "C"
{
#include <blasfeo.h>
}
#include "../AUX/FatropMemory.hpp"
#include "../AUX/FatropLinearAlgebra.hpp"
#include "../AUX/FatropVector.hpp"
#if DEBUG
#include <assert.h>
#endif
using namespace std;
namespace fatrop
{
    // for debuggin purposes
    void fatrop_potrf_l_mn(int m, int n, struct blasfeo_dmat *sC, int ci, int cj, struct blasfeo_dmat *sD, int di, int dj);
    void test();
    /** \brief D <= alpha * B * A^{-1} , with A lower triangular employing explicit inverse of diagonal, fatrop uses its own (naive) implementation since it  not implemented yet in blasfeo*/
    void fatrop_dtrsm_rlnn(int m, int n, double alpha, MAT *sA, int offs_ai, int offs_aj, MAT *sB, int offs_bi, int offs_bj, MAT *sD, int offs_di, int offs_dj);
    // /** \brief D <= alpha * B * A^{-1} , with A lower triangular employing explicit inverse of diagonal, fatrop uses its own (naive) implementation since it  not implemented yet in blasfeo*/
    // void fatrop_dtrsm_rlnn_alt(int m, int n, double alpha, MAT *sA, int offs_ai, int offs_aj, MAT *sB, int offs_bi, int offs_bj, MAT *sD, int offs_di, int offs_dj);
    // B <= B + alpha*A^T (B is mxn)
    void fatrop_dgead_transposed(int m, int n, double alpha, struct blasfeo_dmat *sA, int offs_ai, int offs_aj, struct blasfeo_dmat *sB, int offs_bi, int offs_bj);
    /** this class is used for blasfeo matrices*/
    class fatrop_matrix_bf : public fatrop_matrix
    {
    public:
        /** \brief constructor memory still has to be allocated*/
        fatrop_matrix_bf(const int nrows, const int ncols, const int row_offset, const int col_offset) : row_offset_(row_offset), col_offset_(col_offset), nrows_(nrows), ncols_(ncols) {}
        /** \brief constructor memory already allocated*/
        fatrop_matrix_bf(const int nrows, const int ncols, const int row_offset, const int col_offset, MAT *matbf) : mat_(matbf), row_offset_(row_offset), col_offset_(col_offset), nrows_(nrows), ncols_(ncols) {}
        /** \brief type conversion to blasfeo matrix pointer*/
        inline explicit operator MAT *() const
        {
            return this->mat_;
        }
        /** \brief acces to element of matrix */
        double &at(const int ai, const int aj) const
        {
#if DEBUG
            assert(ai < nrows_);
            assert(aj < ncols_);
#endif
            return MATEL(mat_, ai + row_offset_, aj + col_offset_);
        };
        /** \brief get element of matrix */
        double get_el(const int ai, const int aj) const { return this->at(ai, aj); };
        /** \brief get number of rows */
        int nrows() const { return nrows_; };
        /** \brief get number of cols */
        int ncols() const { return ncols_; };
        /** \brief copies all elements from a given fatrop_matrix to this matrix*/
        void operator=(const fatrop_matrix &fm)
        {
            for (int ai = 0; ai < fm.nrows(); ai++)
            // for (int ai = 0; ai < nrows_; ai++)
            {
                for (int aj = 0; aj < fm.ncols(); aj++)
                // for (int aj = 0; aj < ncols_; aj++)
                {
                    this->at(ai, aj) = fm.get_el(ai, aj);
                }
            }
        }
        /** \brief set data pointer*/
        void set_datap(MAT *matbf)
        {
            mat_ = matbf;
        }
        /** \brief take a block of size (p,q), starting at (i,j)*/
        fatrop_matrix_bf block(const int i, const int j, const int p, const int q) const
        {
            return fatrop_matrix_bf(p, q, row_offset_ + i, col_offset_ + j, this->mat_);
        }

    private:
        MAT *mat_ = NULL;
        const int row_offset_;
        const int col_offset_;
        const int nrows_;
        const int ncols_;
    };

    /** \brief this class is used for the allocation of a blasfeo matrix, the dimsensions are set from a vector */
    class fatrop_memory_matrix_bf : public fatrop_memory_el_base
    {
    public:
        /** \brief constuction for allocation on fatrop_memory_allocator*/
        fatrop_memory_matrix_bf(const FatropVector<int> &nrows, const FatropVector<int> &ncols, int N, fatrop_memory_allocator &fma) : N_(N), nrows_(nrows), ncols_(ncols)
        // TODO: if rvalue-reference is used -> unecessary copy, use move sementics instead.
        {
            fma.add(*this);
        }
        fatrop_memory_matrix_bf(const int nrows, const int ncols, int N, fatrop_memory_allocator &fma) : N_(N), nrows_(vector<int>(N, nrows)), ncols_(vector<int>(N, ncols))
        {
            fma.add(*this);
        }
        /** \brief calculate memory size*/
        int memory_size() const
        {
            int result = 0;
            // size to store structs
            result += N_ * sizeof(MAT);
            // sufficient space for cache alignment
            result = (result + LEVEL1_DCACHE_LINE_SIZE - 1) / LEVEL1_DCACHE_LINE_SIZE * LEVEL1_DCACHE_LINE_SIZE + LEVEL1_DCACHE_LINE_SIZE;
            // size to store date
            for (int i = 0; i < N_; i++)
            {
                result += MEMSIZE_MAT(nrows_.at(i), ncols_.at(i));
            }
            return result;
        };
        /** \brief set up memory element and advance pointer */
        void set_up(char *&data_p)
        {
            MAT *bf_ptr = (MAT *)data_p;
            this->mat = bf_ptr;
            bf_ptr += N_;
            // align with cache line
            long long l_ptr = (long long)bf_ptr;
            l_ptr = (l_ptr + LEVEL1_DCACHE_LINE_SIZE - 1) / LEVEL1_DCACHE_LINE_SIZE * LEVEL1_DCACHE_LINE_SIZE;
            data_p = (char *)l_ptr;
            double *d_ptr_begin = (double *)data_p;
            for (int i = 0; i < N_; i++)
            {
                CREATE_MAT(nrows_.at(i), ncols_.at(i), mat + i, data_p);
                data_p += (mat + i)->memsize;
            }
            double *d_ptr_end = (double *)data_p;
            for (double *d_ptr = d_ptr_begin; d_ptr < d_ptr_end; d_ptr++)
            {
                *d_ptr = 0.0;
            }
        }
        /** \brief get fatrop matrix bf */
        fatrop_matrix_bf operator[](const int N) const
        {
#if DEBUG
            assert(N < N_);
#endif
            MAT *resmat = mat + N;
            fatrop_matrix_bf res(resmat->m, resmat->n, 0, 0, resmat);
            return res;
        }
        /** \brief get first blasfeo_xmat* struct */
        explicit operator MAT *() const
        {
            return mat;
        }

    private:
        MAT *mat;
        const int N_;
        const FatropVector<int> nrows_;
        const FatropVector<int> ncols_;
    };
    /** this class is used for blasfeo vectors*/
    class fatrop_vector_bf : public fatrop_vector
    {
    public:
        /** \brief constructor memory still has to be allocated*/
        fatrop_vector_bf(const int nels, const int offset) : offset_(offset), nels_(nels) {}
        /** \brief constructor memory already allocated*/
        fatrop_vector_bf(const int nels, const int offset, VEC *vecbf) : vec_(vecbf), offset_(offset), nels_(nels) {}
        /** \brief type conversion to blasfeo vector pointer*/
        inline explicit operator VEC *() const
        {
            return this->vec_;
        }
        /** \brief acces to element of matrix */
        double &at(const int ai) const
        {
#if DEBUG
            assert(ai < nels_);
#endif
            return VECEL(vec_, ai + offset_);
        };
        /** \brief get element of vector */
        double get_el(const int ai) const { return this->at(ai); };
        /** \brief get number of elements */
        int nels() const { return nels_; };
        /** \brief copies all elements from a given fatrop_vector to this vector*/
        void operator=(const fatrop_vector &fm)
        {
            for (int ai = 0; ai < nels_; ai++)
            {
                this->at(ai) = fm.get_el(ai);
            }
        }
        /** \brief set data pointer*/
        void set_datap(VEC *vecbf)
        {
            vec_ = vecbf;
        }
        /** \brief take a block of size (p), starting at (i)*/
        fatrop_vector_bf block(const int i, const int p) const
        {
            return fatrop_vector_bf(p, offset_ + i, this->vec_);
        }

    private:
        VEC *vec_ = NULL;
        const int offset_;
        const int nels_;
    };
    /** \brief this class is used for the allocation of a blasfeo vector, the dimsensions are set from a vector */
    class fatrop_memory_vector_bf : public fatrop_memory_el_base
    {
    public:
        /** \brief constuction for allocation on fatrop_memory_allocator*/
        fatrop_memory_vector_bf(const FatropVector<int> &nels, int N, fatrop_memory_allocator &fma) : N_(N), nels_(nels)
        // TODO: if rvalue-reference is used -> unecessary copy, use move sementics instead.
        {
            fma.add(*this);
        }
        fatrop_memory_vector_bf(const int nels, int N, fatrop_memory_allocator &fma) : N_(N), nels_(vector<int>(N, nels))
        {
            fma.add(*this);
        }
        /** \brief calculate memory size*/
        int memory_size() const
        {
            int result = 0;
            // size to store structs
            result += N_ * sizeof(VEC);
            // sufficient space for cache alignment
            result = (result + LEVEL1_DCACHE_LINE_SIZE - 1) / LEVEL1_DCACHE_LINE_SIZE * LEVEL1_DCACHE_LINE_SIZE + LEVEL1_DCACHE_LINE_SIZE;
            // size to store date
            for (int i = 0; i < N_; i++)
            {
                result += MEMSIZE_VEC(nels_.at(i));
            }
            return result;
        };
        /** \brief set up memory element and advance pointer */
        void set_up(char *&data_p)
        {
            VEC *bf_ptr = (VEC *)data_p;
            this->vec = bf_ptr;
            bf_ptr += N_;
            // align with cache line
            long long l_ptr = (long long)bf_ptr;
            l_ptr = (l_ptr + LEVEL1_DCACHE_LINE_SIZE - 1) / LEVEL1_DCACHE_LINE_SIZE * LEVEL1_DCACHE_LINE_SIZE;
            data_p = (char *)l_ptr;
            double *d_ptr_begin = (double *)data_p;
            for (int i = 0; i < N_; i++)
            {
                CREATE_VEC(nels_.at(i), vec + i, data_p);
                data_p += (vec + i)->memsize;
            }
            double *d_ptr_end = (double *)data_p;
            for (double *d_ptr = d_ptr_begin; d_ptr < d_ptr_end; d_ptr++)
            {
                *d_ptr = 0.0;
            }
        }
        /** \brief get fatrop matrix bf */
        fatrop_vector_bf operator[](const int N) const
        {
#if DEBUG
            assert(N < N_);
#endif
            VEC *resvec = vec + N;
            fatrop_vector_bf res(resvec->m, 0, resvec);
            return res;
        }
        /** \brief get first blasfeo_xmat* struct */
        explicit operator VEC *() const
        {
            return vec;
        }

    private:
        VEC *vec;
        const int N_;
        const FatropVector<int> nels_;
    };

    /** \brief this class represents a permutation matrix */
    class fatrop_permutation_matrix : public fatrop_matrix
    {
    public:
        /** \brief constructor memory still has to be allocated */
        fatrop_permutation_matrix(const int dim) : dim_(dim){};
        /** \brief constructor memory already allocated */
        fatrop_permutation_matrix(const int dim, int *data) : dim_(dim), data_(data){};
        /** \brief get number of rows */
        int nrows() const { return dim_; };
        /** \brief get number of columns */
        int ncols() const { return dim_; };
        /** \brief get element of matrix represented by this permutation matrix - only used for debugging and testing purposes */
        double get_el(const int ai, const int aj) const
        {
#if DEBUG
            assert(data_ != NULL);
#endif
            int aj_one = data_[ai];
            int row_curr = ai - 1;
            while (row_curr >= 0)
            {
                if (aj_one == data_[row_curr])
                {
                    aj_one = row_curr;
                }
                row_curr--;
            }
            if (aj == aj_one)
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        };
        void print(const int kmax) const
        {
            for (int k = 0; k < kmax; k++)
            {
                cout << k << " <-> " << data_[k] << endl;
            }
        }
        /** \brief set data pointer*/
        void set_datap(int *data)
        {
            data_ = data;
        }
        /** \brief set data point*/
        void set_datap(const int i, const int val)
        {
#if DEBUG
            assert(data_ != NULL);
            assert(i < dim_);
#endif
            data_[i] = val;
        }
        /** \brief apply row permutation*/
        void PM(const int kmax, MAT *M) const
        {
#if DEBUG
            assert(data_ != NULL);
#endif
            ROWPE(kmax, data_, M);
        }
        /** \brief apply vec permutation*/
        void PV(const int kmax, VEC *V, const int offs) const
        {
#if DEBUG
            assert(data_ != NULL);
#endif
            VECPE(kmax, data_, V, offs);
        }
        /** \brief apply vec permutation*/
        void PtV(const int kmax, VEC *V, const int offs) const
        {
#if DEBUG
            assert(data_ != NULL);
#endif
            VECPEI(kmax, data_, V, offs);
        }
        /** \brief apply row permutation on partial matrix*/
        void PM(const int kmax, const int n, MAT *M, const int ai, const int aj) const
        {
#if DEBUG
            assert(data_ != NULL);
#endif
            // invalidate stored inverse diagonal
            M->use_dA = 0;

            int ii;
            for (ii = 0; ii < kmax; ii++)
            {
                if (data_[ii] != ii)
                    ROWSW(n, M, ai + ii, aj, M, ai + data_[ii], aj);
            }
            return;
        }
        /** \brief apply inverse row permutation*/
        void
        PtM(const int kmax, MAT *M) const
        {
#if DEBUG
            assert(data_ != NULL);
#endif
            ROWPEI(kmax, data_, M);
        }
        /** \brief apply inverse col permutation*/
        void MP(const int kmax, MAT *M) const
        {
#if DEBUG
            assert(data_ != NULL);
#endif
            COLPEI(kmax, data_, M);
        }
        /** \brief apply col permutation*/
        void MPt(const int kmax, MAT *M) const
        {
#if DEBUG
            assert(data_ != NULL);
#endif
            COLPE(kmax, data_, M);
        }
        /** int pointer of permutation vector */
        explicit operator int *() { return data_; };

        // private:
        const int dim_;
        int *data_ = NULL;
    };

    /** \brief this class is used for the allocation of a permutation matrix */
    class fatrop_memory_permutation_matrix : public fatrop_memory_el_base, public fatrop_permutation_matrix
    {
    public:
        /** \brief constructor */
        fatrop_memory_permutation_matrix(const int dim, const int N, fatrop_memory_allocator &fma) : fatrop_permutation_matrix(dim), dim_(dim), N_(N)
        {
            fma.add(*this);
        };
        /** \brief calculate needed memory size*/
        int memory_size() const
        {
            int size = 0;
            size += N_ * sizeof(fatrop_permutation_matrix) + N_ * dim_ * sizeof(int);
            return size;
        }
        /** \brief set up memory*/
        void set_up(char *&char_p)
        {
            perm_p = (fatrop_permutation_matrix *)char_p;
            fatrop_permutation_matrix *perm_pp = perm_p;
            for (int i = 0; i < N_; i++)
            {
                new (perm_pp) fatrop_permutation_matrix(dim_);
                perm_pp++;
            }
            int *data_p = (int *)perm_pp;
            this->set_datap(data_p);
            for (int i = 0; i < N_; i++)
            {
                perm_p[i].set_datap(data_p);
                data_p += dim_;
            }
            char_p = (char *)data_p;
        }
        explicit operator fatrop_permutation_matrix *() { return perm_p; };

    private:
        const int dim_;
        const int N_;
        fatrop_permutation_matrix *perm_p;
    };
    matrix_ind max_el(int m, int n, MAT *matr, int ai, int aj);
#define MAX(a, b)                   \
    (                               \
        {                           \
            __typeof__(a) _a = (a); \
            __typeof__(b) _b = (b); \
            _a > _b ? _a : _b;      \
        })
#define MIN(a, b)                   \
    (                               \
        {                           \
            __typeof__(a) _a = (a); \
            __typeof__(b) _b = (b); \
            _a < _b ? _a : _b;      \
        })
    /** \brief Function to calculate LU factorization result is saved in A, L is unit diagonal */
    void LU_FACT(const int m, const int n, const int n_max, int &rank, MAT *A, PMAT *Pl_p, PMAT *Pr_p, double tol = 1e-12);
    /** \brief Function to calculate LU factorization but A, and result (L and U) are transposed, all indices refer to the dimensions of the original A matrix (and not the transposed one) */
    void LU_FACT_transposed(const int m, const int n, const int n_max, int &rank, MAT *At, PMAT *Pl_p, PMAT *Pr_p, double tol = 1e-12);
    void fatrop_dtrsv_unu(int m, blasfeo_dmat *sA, int ai, int aj, blasfeo_dvec *sx, int xi, blasfeo_dvec *sz, int zi);

} // namespace fatrop
#endif //FATROP_BLASFEO_INCLUDED