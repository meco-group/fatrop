#ifndef FATROP_BLASFEO_INCLUDED
#define FATROP_BLASFEO_INCLUDED

// macros
#define MAT blasfeo_dmat
#define MEMSIZE_MAT blasfeo_memsize_dmat
#define MEMSIZE_MAT blasfeo_memsize_dmat
#define CREATE_MAT blasfeo_create_dmat
#define ROWPE blasfeo_drowpe
#define ROWPEI blasfeo_drowpei
#define COLPE blasfeo_dcolpe
#define COLPEI blasfeo_dcolpei
#define CREATE_MAT blasfeo_create_dmat
#define MATEL BLASFEO_DMATEL

#include <iostream>
extern "C"
{
#include <blasfeo.h>
}
#include "FatropMemory.hpp"
#include "FatropLinearAlgebra.hpp"
#if DEBUG
#include <assert.h>
#endif
using namespace std;
namespace fatrop
{
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
            for (int ai = 0; ai < nrows_; ai++)
            {
                for (int aj = 0; aj < ncols_; aj++)
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
    /** \brief this class is used for the allocation of a blasfeo matrix */
    class fatrop_memory_matrix_bf : public fatrop_memory_el_base, public fatrop_matrix_bf
    {
    public:
        /** \brief constuction for allocation on fatrop_memory_allocator*/
        fatrop_memory_matrix_bf(int nrows, int ncols, int N, fatrop_memory_allocator &fma) : fatrop_matrix_bf(nrows, ncols, 0, 0), N_(N), nrows_(nrows), ncols_(ncols)
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
            result += N_ * MEMSIZE_MAT(nrows_, ncols_);
            return result;
        };
        /** \brief set up memory element and advance pointer */
        void set_up(char *&data_p)
        {
            MAT *bf_ptr = (MAT *)data_p;
            this->mat = bf_ptr;
            fatrop_matrix_bf::set_datap(bf_ptr);
            bf_ptr += N_;
            // align with cache line
            long long l_ptr = (long long)bf_ptr;
            l_ptr = (l_ptr + LEVEL1_DCACHE_LINE_SIZE - 1) / LEVEL1_DCACHE_LINE_SIZE * LEVEL1_DCACHE_LINE_SIZE;
            data_p = (char *)l_ptr;
            double *d_ptr_begin = (double *)data_p;
            for (int i = 0; i < N_; i++)
            {
                CREATE_MAT(nrows_, ncols_, mat + i, data_p);
                data_p += (mat + i)->memsize;
            }
            double *d_ptr_end = (double *)data_p;
            for (double *d_ptr = d_ptr_begin; d_ptr < d_ptr_end; d_ptr++)
            {
                *d_ptr = 0.0;
            }
        }
        /** \brief copy operator */
        void operator=(const fatrop_matrix &fm)
        {
            fatrop_matrix_bf::operator=(fm);
        }

    private:
        MAT *mat;
        const int N_;
        const int nrows_;
        const int ncols_;
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
        /** \brief apply inverse row permutation*/
        void PtM(const int kmax, MAT *M) const
        {
#if DEBUG
            assert(data_ != NULL);
#endif
            ROWPEI(kmax, data_, M);
        }
        /** \brief apply col permutation*/
        void MP(const int kmax, MAT *M) const
        {
#if DEBUG
            assert(data_ != NULL);
#endif
            COLPE(kmax, data_, M);
        }
        /** \brief apply inverse col permutation*/
        void MPt(const int kmax, MAT *M) const
        {
#if DEBUG
            assert(data_ != NULL);
#endif
            COLPEI(kmax, data_, M);
        }

    private:
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
            size += N_ * dim_ * sizeof(int);
            return size;
        }
        /** \brief set up memory*/
        void set_up(char *&char_p)
        {
            data_ = (int *)char_p;
            this->set_datap(data_);
            char_p += memory_size();
        }
        /** \brief get n-th permutation vector pointer*/
        int *perm_vector(const int n) const
        {
#if DEBUG
            assert(n < N_);
#endif
            return data_ + dim_ * n;
        }

    private:
        const int dim_;
        const int N_;
        int *data_;
    };
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
    void LU_FACT(const int m, const int n, const int n_max, int &rank, MAT *A, const int *perm_left, const int *perm_right){
        int minmn = MIN(m,n_max);
        int j =0;
        for(int i=0; i<minmn; i++){
            matrix_ind max_curr = max_el(m,n_max, A,i,i);
        }
    };
} // namespace fatrop
#endif //FATROP_BLASFEO_INCLUDED