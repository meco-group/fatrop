
#ifndef FATROP_BLASFEO_INCLUDED
#define FATROP_BLASFEO_INCLUDED

// macros
#define MAT blasfeo_dmat
#define MEMSIZE_MAT blasfeo_memsize_dmat
#define MEMSIZE_MAT blasfeo_memsize_dmat
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
    class fatrop_matrix_bf : public fatrop_matrix
    {
    public:
        fatrop_matrix_bf(const int nrows, const int ncols, const int row_offset, const int col_offset) : row_offset_(row_offset), col_offset_(col_offset), nrows_(nrows), ncols_(ncols) {}
        fatrop_matrix_bf(const int nrows, const int ncols, const int row_offset, const int col_offset, MAT *matbf) : mat_(matbf), row_offset_(row_offset), col_offset_(col_offset), nrows_(nrows), ncols_(ncols){}
        inline explicit operator MAT *() const
        {
            return this->mat_;
        }
        double get_el(const int ai, const int aj) const { return MATEL(mat_, ai, aj); };
        /** \brief acces to element of matrix */
        double &at(const int ai, const int aj) const { return MATEL(mat_, ai, aj); };
        int nrows() const { return nrows_; };
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
        void set_datap(MAT *matbf)
        {
            mat_ = matbf;
        }

    private:
        MAT *mat_ = NULL;
        const int row_offset_;
        const int col_offset_;
        const int nrows_;
        const int ncols_;
    };
    class fatrop_memory_matrix : public fatrop_memory_el_base, public fatrop_matrix_bf
    {
    public:
        /** \brief constuction for allocation on fatrop_memory_allocator*/
        fatrop_memory_matrix(int nrows, int ncols, int N, fatrop_memory_allocator &fma) : fatrop_matrix_bf(nrows, ncols, 0, 0), N_(N), nrows_(nrows), ncols_(ncols)
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
    class fatrop_permutation_matrix : public fatrop_matrix
    {
    public:
        fatrop_permutation_matrix(const int dim) : dim_(dim){};
        fatrop_permutation_matrix(const int dim, int *data) : dim_(dim), data_(data){};
        int nrows() const { return dim_; };
        int ncols() const { return dim_; };
        double get_el(const int ai, const int aj) const
        {
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
        void set_datap(int *data)
        {
            data_ = data;
        }

    private:
        const int dim_;
        int *data_ = NULL;
    };

    /** \brief this class represents a permutation matrix */
    class fatrop_memory_permutation_matrix : public fatrop_memory_el_base, public fatrop_permutation_matrix
    {
    public:
        fatrop_memory_permutation_matrix(const int dim, const int N, fatrop_memory_allocator &fma) : fatrop_permutation_matrix(dim), dim_(dim), N_(N)
        {
            fma.add(*this);
        };
        int memory_size() const
        {
            int size = 0;
            size += N_ * dim_ * sizeof(int);
            return size;
        }
        void set_up(char *&char_p)
        {
            data_ = (int *)char_p;
            this->set_datap(data_);
            char_p += memory_size();
        }
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
} // namespace fatrop
#endif //FATROP_BLASFEO_INCLUDED