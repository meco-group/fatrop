//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

/**
 * @file vector_index.hppp
 * @brief Defines the IndexVector class and related functionality for efficient index manipulation.
 */

#ifndef __fatrop_common_index_vector__
#define __fatrop_common_index_vector__
#include <cstddef>
#include <numeric>
#include <span>
#include <vector>

/**
 * @namespace fatrop
 */
namespace fatrop
{
    template <typename Derived, typename Index> class VecIndexBlock;

    template <typename Derived, typename Index> class VecIndex
    {
    public:
        /**
         * @brief Accesses the element at the given index.
         *
         * @param i Index of the element.
         * @return Index The value at index i.
         */
        Index operator[](const Index i) const
        {
            return static_cast<const Derived *>(this)->operator[](i);
        }

        /**
         * @brief Gets the size (number of elements) of the vector.
         *
         * @return Index The size of the vector.
         */
        Index m() const { return static_cast<const Derived *>(this)->m(); }

        /**
         * @brief Extracts a sub-block (segment) of the vector.
         *
         * @param size Size of the block.
         * @param start Starting index of the block.
         * @return VecIndexBlock<Derived> The sub-block of the vector.
         */
        VecIndexBlock<Derived, Index> block(Index start, Index size) const
        {
            return VecIndexBlock<Derived, Index>(*static_cast<const Derived *>(this), start, size);
        }

        // Various mathematical operations defined as friend functions
        // They allow mathematical operations like sum, norms, and transformations on
        // vectors.
        friend Index sum(const VecIndex<Derived, Index> &vec)
        {
            Index ret = 0;
            for (Index i = 0; i < vec.m(); i++)
            {
                ret += vec[i];
            }
            return ret;
        }
    };

    template <typename Derived, typename Index>
    class VecIndexBlock : public VecIndex<VecIndexBlock<Derived, Index>, Index>
    {
    public:
        VecIndexBlock(const Derived &a, const Index start, const Index size) : a(a), ai_(start), m_(size) {};
        Index operator[](const Index i) const { return a[ai_ + i]; }
        Index m() const { return m_; }
    private:
        const Derived &a;
        const Index ai_;
        const Index m_;
    };

    template <typename Index>
    class VecIndexAllocated;

    template <typename Index>
    class IndexVecView: public VecIndex<IndexVecView<Index>, Index>
    {
    public:
        IndexVecView(VecIndexAllocated<Index>& vec, const Index start, const Index size):  vec_(vec), ai_(start), m_(size) {};
        inline Index operator[](const Index i) const;
        inline Index m() const;
        IndexVecView<Index> block(Index start, Index size) const
        {
            return IndexVecView<Index>(vec_, ai_ + start, size);
        }
        template <typename Derived>
        IndexVecView<Index>& operator=(const VecIndex<Derived, Index> &vec_in)
        {
            fatrop_dbg_assert(this->size() == vec_in.m() &&
                              "Vectors must be same size for asignment");
            for (Index i = 0; i < m(); i++)
            {
                this->operator[](i) = vec_in(i);
            }
            return *this;
        }
        IndexVecView<Index>& operator=(const IndexVecView& vec_in)
        {
            *this = static_cast<const VecIndex<IndexVecView<Index>, Index>&>(vec_in);
        }
        VecIndexAllocated<Index>& vec_;
        Index ai_;
        Index m_;
    };

    template <typename Index>
    class VecIndexAllocated : public IndexVecView<Index>,
                              public std::vector<Index>
    {
    public:
        // using std::vector<Index>::vector;
        VecIndexAllocated(std::vector<Index>&& vec): IndexVecView<Index>(*this, 0, vec.size()), std::vector<Index>(std::move(vec)) {};
        VecIndexAllocated(const std::vector<Index>& vec): IndexVecView<Index>(*this, 0, vec.size()), std::vector<Index>(vec) {};
        using std::vector<Index>::operator[];
        template <typename Derived>
        Index m() const { return this->size(); }
    };
    // implementation

    template <typename Index>
        Index IndexVecView<Index>::operator[](const Index i) const { return vec_[i+ai_]; }
    template <typename Index>
        Index IndexVecView<Index>::m() const { return m_; }
} // namespace fatrop

#endif //__fatrop_common_index_vector__
