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
#define GESC blasfeo_dgesc
#define VECCP blasfeo_dveccp
#define VECCPSC blasfeo_dveccpsc
#define VECCPR fatrop_dveccp_reversed
#define TRSM_RLNN fatrop_dtrsm_rlnn_alt // TODO this is not implemented by blasfeo so we defined our own (naive) implementation
// #define TRSM_LLNN blasfeo_dtrsm_llnn
#define TRSM_RLTN blasfeo_dtrsm_rltn
#define VECEL BLASFEO_DVECEL
#define MEMSIZE_VEC blasfeo_memsize_dvec
#define CREATE_VEC blasfeo_create_dvec
#define GEMM_NT blasfeo_dgemm_nt
#define GEAD blasfeo_dgead
#define SYRK_LN_MN blasfeo_dsyrk_ln_mn
#define GETR blasfeo_dgetr
#define TRTR_L blasfeo_dtrtr_l
#define POTRF_L_MN blasfeo_dpotrf_l_mn
#define ROWEX blasfeo_drowex
#define ROWIN blasfeo_drowin
#define ROWAD fatrop_drowad
#define TRSV_LTN blasfeo_dtrsv_ltn
#define TRSV_LNN blasfeo_dtrsv_lnn
#define TRSV_UNU fatrop_dtrsv_unu
#define TRSV_UTN blasfeo_dtrsv_utn
#define TRSV_UTU fatrop_dtrsv_utu
#define GEMV_T blasfeo_dgemv_t
#define GEMV_N blasfeo_dgemv_n
#define VECSE blasfeo_dvecse
#define VECSC blasfeo_dvecsc
#define PACKMAT blasfeo_pack_dmat
#define UNPACKVEC blasfeo_unpack_dvec
#define PACKVEC blasfeo_pack_dvec
#define PMAT PermMat
#define AXPY blasfeo_daxpy
#define AXPBY blasfeo_daxpby
#define DOT blasfeo_ddot
#define GESE blasfeo_dgese
#define DIARE blasfeo_ddiare
#define COLSC blasfeo_dcolsc
#define VECMUL blasfeo_dvecmul
#define VECMULACC blasfeo_dvecmulacc

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
#include <iostream>
extern "C"
{
#include <blasfeo.h>
}
#include "auxiliary/LinearAlgebra.hpp"
#include "auxiliary/FatropVector.hpp"
#include "auxiliary/Common.hpp"
#if DEBUG
#include <assert.h>
#endif
namespace fatrop
{
    void fatrop_dcolsc(int kmax, double alpha, struct blasfeo_dmat *sA, int ai, int aj);
    // copy elements from sx to sy but in reversed order to avoid aliasing issues in recursion
    void fatrop_dveccp_reversed(int m, struct blasfeo_dvec *sx, int xi, struct blasfeo_dvec *sy, int yi);
    // for debugging purposes
    // void fatrop_potrf_l_mn(int m, int n, struct blasfeo_dmat *sC, int ci, int cj, struct blasfeo_dmat *sD, int di, int dj);
    void test();
    /** \brief D <= alpha * B * A^{-1} , with A lower triangular employing explicit inverse of diagonal, fatrop uses its own (naive) implementation since it  not implemented yet in blasfeo */
    void fatrop_dtrsm_rlnn(int m, int n, double alpha, MAT *sA, int offs_ai, int offs_aj, MAT *sB, int offs_bi, int offs_bj, MAT *sD, int offs_di, int offs_dj);
    /** \brief D <= alpha * B * A^{-1} , with A lower triangular employing explicit inverse of diagonal, fatrop uses its own (naive) implementation since it  not implemented yet in blasfeo */
    void fatrop_dtrsm_rlnn_alt(int m, int n, double alpha, MAT *sA, int offs_ai, int offs_aj, MAT *sB, int offs_bi, int offs_bj, MAT *sD, int offs_di, int offs_dj);
    /** \brief B <= B + alpha*A^T (B is mxn) */
    void fatrop_dgead_transposed(int m, int n, double alpha, struct blasfeo_dmat *sA, int offs_ai, int offs_aj, struct blasfeo_dmat *sB, int offs_bi, int offs_bj);
    void fatrop_identity(const int m, MAT *sA, const int ai, const int aj);
    void fatrop_drowad(int kmax, double alpha, struct blasfeo_dvec *sx, int xi, struct blasfeo_dmat *sA, int ai, int aj);
    /** \brief this class is used for blasfeo matrices*/
    class FatropMatBF : public FatropMat
    {
    public:
        /** \brief constructor memory still has to be allocated*/
        FatropMatBF(const int nrows, const int ncols, const int row_offset, const int col_offset);
        /** \brief constructor memory already allocated*/
        FatropMatBF(const int nrows, const int ncols, const int row_offset, const int col_offset, MAT *matbf);
        /** \brief constructor memory already allocated*/
        FatropMatBF(MAT *matbf);
        /** \brief type conversion to blasfeo matrix pointer*/
        inline explicit operator MAT *() const
        {
            return this->mat_;
        }
        /** \brief acces to element of matrix */
        inline double &at(const int ai, const int aj) const
        {
#if DEBUG
            assert(ai < nrows_);
            assert(aj < ncols_);
#endif
            return MATEL(mat_, ai + row_offset_, aj + col_offset_);
        };
        /** \brief get element of matrix */
        inline double get_el(const int ai, const int aj) const { return this->at(ai, aj); };
        /** \brief get number of rows */
        inline int nrows() const { return nrows_; };
        /** \brief get number of cols */
        inline int ncols() const { return ncols_; };
        /** \brief copies all elements from a given fatrop_matrix to this matrix*/
        void operator=(const FatropMat &fm);
        /** \brief set data pointer*/
        void set_datap(MAT *matbf)
        {
            mat_ = matbf;
        }
        /** \brief take a block of size (p,q), starting at (i,j)*/
        FatropMatBF block(const int i, const int j, const int p, const int q) const
        {
            return FatropMatBF(p, q, row_offset_ + i, col_offset_ + j, this->mat_);
        }

    private:
        MAT *mat_ = NULL;
        const int row_offset_;
        const int col_offset_;
        const int nrows_;
        const int ncols_;
    };

    /** \brief this class is used for the allocation of a blasfeo matrix, the dimsensions are set from a vector */
    class FatropMemoryMatBF
    {
    public:
        /** \brief constuction for allocation on fatrop_memory_allocator*/
        FatropMemoryMatBF(const FatropVector<int> &nrows, const FatropVector<int> &ncols, int N);
        // TODO: if rvalue-reference is used -> unecessary copy, use move sementics instead.;
        FatropMemoryMatBF(const int nrows, const int ncols, int N);
        /** \brief calculate memory size*/
        int memory_size() const;
        /** \brief set up memory element and advance pointer */
        void set_up();
        /** \brief get fatrop matrix bf */
        FatropMatBF operator[](const int N) const;
        /** \brief get first blasfeo_xmat* struct */
        explicit operator MAT *() const
        {
            return mat;
        }
        FatropMemoryMatBF(const FatropMemoryMatBF &cpy) = delete;
        FatropMemoryMatBF &operator=(const FatropMemoryMatBF &) = delete;
        ~FatropMemoryMatBF();

    private:
        void *mem = NULL;
        MAT *mat;
        const int N_;
        const FatropVector<int> nrows_;
        const FatropVector<int> ncols_;
    };
    /** this class is used for blasfeo vectors*/
    class FatropVecBF : public FatropVec
    {
    public:
        /** \brief constructor memory still has to be allocated*/
        FatropVecBF(const int nels, const int offset);
        /** \brief constructor memory already allocated*/
        FatropVecBF(const int nels, const int offset, VEC *vecbf);
        /** \brief type conversion to blasfeo vector pointer*/
        explicit operator VEC *() const;
        /** \brief access to element of matrix */
        double &at(const int ai) const;
        /** \brief get element of vector */
        double get_el(const int ai) const;
        /** \brief get number of elements */
        int nels() const;
        /** \brief get offset */
        int offset() const;
        /** \brief copies all elements from a given fatrop_vector to this vector*/
        void operator=(const FatropVec &fm);
        void copy(const FatropVecBF &fm);
        void copyto(std::vector<double>& dest) const;
        void operator=(const std::vector<double> &fm);
        /** \brief set data pointer*/
        void set_datap(VEC *vecbf);
        /** \brief take a block of size (p), starting at (i)*/
        FatropVecBF block(const int i, const int p) const;
        void SwapWith(FatropVecBF &vb);
        void SetConstant(double constant) const;

    protected:
        VEC *vec_ = NULL;
        const int offset_;
        const int nels_;
    };

    void axpy(const double alpha, const FatropVecBF &va, const FatropVecBF &vb, const FatropVecBF &vc);
    void copy(const FatropVecBF &va, const FatropVecBF &vb);
    void axpby(const double alpha, const FatropVecBF &va, const double beta, const FatropVecBF &vb, const FatropVecBF &vc);
    double dot(const FatropVecBF &va, FatropVecBF &vb);
    double Linf(const FatropVecBF &va);
    double LinfScaled(const FatropVecBF &va, const FatropVecBF &scales);
    double minabs(const FatropVecBF &va);
    double L1(const FatropVecBF &va);

    /** \brief this class is used for the allocation of a blasfeo vector, the dimsensions are set from a vector */
    class FatropMemoryVecBF
    {
    public:
        /** \brief constuction for allocation on MemoryAllocator*/
        FatropMemoryVecBF(const FatropVector<int> &nels, int N);
        // TODO: if rvalue-reference is used -> unecessary copy, use move sementics instead.;
        FatropMemoryVecBF(const int nels, int N);
        /** \brief calculate memory size*/
        int memory_size() const;
        /** \brief set up memory element and advance pointer */
        void set_up();
        /** \brief get fatrop matrix bf */
        FatropVecBF operator[](const int N) const;
        /** \brief get first blasfeo_xmat* struct */
        explicit operator VEC *() const
        {
            return vec;
        }
        FatropMemoryVecBF(const FatropMemoryVecBF &cpy) = delete;
        FatropMemoryVecBF &operator=(const FatropMemoryVecBF &) = delete;
        ~FatropMemoryVecBF();

    private:
        void *mem = NULL;
        VEC *vec;
        const int N_;
        const FatropVector<int> nels_;
    };

    /** \brief this class represents a permutation matrix */
    class PermMat : public FatropMat
    {
    public:
        /** \brief constructor memory still has to be allocated */
        PermMat(const int dim);
        ;
        /** \brief constructor memory already allocated */
        PermMat(const int dim, int *data);
        ;
        /** \brief get number of rows */
        int nrows() const { return dim_; };
        /** \brief get number of columns */
        int ncols() const { return dim_; };
        /** \brief get element of matrix represented by this permutation matrix - only used for debugging and testing purposes */
        double get_el(const int ai, const int aj) const;
        void print(const int kmax) const;
        /** \brief set data pointer*/
        void set_datap(int *data);
        /** \brief set data point*/
        void set_datap(const int i, const int val);
        /** \brief apply row permutation*/
        void PM(const int kmax, MAT *M) const;
        /** \brief apply vec permutation*/
        void PV(const int kmax, VEC *V, const int offs) const;
        /** \brief apply vec permutation*/
        void PtV(const int kmax, VEC *V, const int offs) const;
        /** \brief apply row permutation on partial matrix*/
        void PM(const int kmax, const int n, MAT *M, const int ai, const int aj) const;
        /** \brief apply inverse row permutation*/
        void
        PtM(const int kmax, MAT *M) const;
        /** \brief apply inverse col permutation*/
        void MP(const int kmax, MAT *M) const;
        /** \brief apply col permutation*/
        void MPt(const int kmax, MAT *M) const;
        /** int pointer of permutation vector */
        explicit operator int *() { return data_; };

        // private:
        const int dim_;
        int *data_ = NULL;
    };

    /** \brief this class is used for the allocation of a permutation matrix */
    class MemoryPermMat : public PermMat
    {
    public:
        /** \brief constructor */
        MemoryPermMat(const int dim, const int N);
        /** \brief calculate needed memory size*/
        int memory_size() const;
        /** \brief set up memory*/
        void set_up();
        explicit operator PermMat *()
        {
            return perm_p;
        };
        MemoryPermMat(const MemoryPermMat &cpy) = delete;
        MemoryPermMat &operator=(const MemoryPermMat &) = delete;
        ~MemoryPermMat();

    private:
        void *mem = NULL;
        const int dim_;
        const int N_;
        PermMat *perm_p;
    };
    MatrixInd max_el(int m, int n, MAT *matr, int ai, int aj);

    /** \brief Function to calculate LU factorization result is saved in A, L is lower unitriangular */
    void LU_FACT(const int m, const int n, const int n_max, int &rank, MAT *A, PMAT *Pl_p, PMAT *Pr_p, double tol = 1e-8);
    /** \brief Function to calculate LU factorization but A, and result (L and U) are transposed, all indices refer to the dimensions of the original A matrix (and not the transposed one) */
    void LU_FACT_transposed(const int m, const int n, const int n_max, int &rank, MAT *At, PMAT *Pl_p, PMAT *Pr_p, double tol = 1e-5);
    void fatrop_dtrsv_unu(const int m, const int n, blasfeo_dmat *sA, const int ai, const int aj, blasfeo_dvec *sx, const int xi, blasfeo_dvec *sz, const int zi);
    void fatrop_dtrsv_utu(const int m, blasfeo_dmat *sA, const int ai, const int aj, blasfeo_dvec *sx, const int xi, blasfeo_dvec *sz, const int zi);
} // namespace fatrop
#endif // FATROP_BLASFEO_INCLUDED