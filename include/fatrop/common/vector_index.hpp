//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

/**
 * @file vector_index.hpp
 * @brief Defines classes and functionality for efficient index vector manipulation.
 * 
 * This file contains template classes for index vectors, including VecIndex, VecIndexBlock,
 * IndexVecView, and VecIndexAllocated. These classes provide a flexible and efficient way
 * to work with vectors of indices, supporting operations like sub-block extraction,
 * element access, and size queries.
 */

#ifndef __fatrop_common_index_vector__
#define __fatrop_common_index_vector__
#include <cstddef>
#include <numeric>
#include <span>
#include <vector>
#include "fatrop/common/exception.hpp"

/**
 * @namespace fatrop
 */
namespace fatrop
{
    /**
     * @brief Forward declaration of VecIndexBlock class.
     * 
     * @tparam Derived The derived class type for CRTP.
     * @tparam Index The index type used in the vector.
     */
    template <typename Derived, typename Index> class VecIndexBlock;

    /**
     * @brief Base class for index vectors using CRTP (Curiously Recurring Template Pattern).
     * 
     * This class provides a common interface for different types of index vectors,
     * allowing for polymorphic behavior without virtual functions.
     * 
     * @tparam Derived The derived class type.
     * @tparam Index The index type used in the vector.
     */
    template <typename Derived, typename Index> class VecIndex
    {
    public:
        /**
         * @brief Accesses the element at the given index.
         *
         * This operator provides read-only access to the elements of the vector.
         *
         * @param i Index of the element to access.
         * @return Index The value at index i.
         */
        Index operator[](const Index i) const
        {
            return static_cast<const Derived *>(this)->operator[](i);
        }

        /**
         * @brief Gets the size (number of elements) of the vector.
         *
         * This method returns the total number of elements in the vector.
         *
         * @return Index The size of the vector.
         */
        Index m() const { return static_cast<const Derived *>(this)->m(); }

        /**
         * @brief Extracts a sub-block (segment) of the vector.
         *
         * This method allows for creating a view of a portion of the vector
         * without copying the data.
         *
         * @param start Starting index of the block.
         * @param size Size of the block.
         * @return VecIndexBlock<Derived, Index> A view representing the sub-block of the vector.
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

    /**
     * @brief Represents a block (sub-vector) of an index vector.
     * 
     * This class provides a view into a portion of a larger vector without copying data.
     * 
     * @tparam Derived The type of the original vector.
     * @tparam Index The index type used in the vector.
     */
    template <typename Derived, typename Index>
    class VecIndexBlock : public VecIndex<VecIndexBlock<Derived, Index>, Index>
    {
    public:
        /**
         * @brief Constructs a VecIndexBlock.
         * 
         * @param a The original vector.
         * @param start The starting index of the block in the original vector.
         * @param size The size of the block.
         */
        VecIndexBlock(const Derived &a, const Index start, const Index size) : a(a), ai_(start), m_(size) {};
        
        /**
         * @brief Accesses an element in the block.
         * 
         * @param i Index within the block.
         * @return Index The value at the specified index.
         */
        Index operator[](const Index i) const { return a[ai_ + i]; }
        
        /**
         * @brief Returns the size of the block.
         * 
         * @return Index The number of elements in the block.
         */
        Index m() const { return m_; }
    private:
        const Derived &a;  ///< Reference to the original vector
        const Index ai_;   ///< Starting index of the block in the original vector
        const Index m_;    ///< Size of the block
    };

    template <typename Index>
    class VecIndexAllocated;

    /**
     * @brief Provides a view into a VecIndexAllocated object.
     * 
     * This class allows for creating a view of a portion of a VecIndexAllocated
     * without copying the data.
     * 
     * @tparam Index The index type used in the vector.
     */
    template <typename Index>
    class IndexVecView: public VecIndex<IndexVecView<Index>, Index>
    {
    public:
        /**
         * @brief Constructs an IndexVecView.
         * 
         * @param vec Reference to the VecIndexAllocated object.
         * @param start Starting index of the view.
         * @param size Size of the view.
         */
        IndexVecView(VecIndexAllocated<Index>& vec, const Index start, const Index size):  vec_(vec), ai_(start), m_(size) {};
        
        inline Index operator[](const Index i) const;
        inline Index m() const;
        
        /**
         * @brief Creates a sub-block of this view.
         * 
         * @param start Starting index within this view.
         * @param size Size of the new block.
         * @return IndexVecView<Index> A new view representing the sub-block.
         */
        IndexVecView<Index> block(Index start, Index size) const
        {
            return IndexVecView<Index>(vec_, ai_ + start, size);
        }
        
        /**
         * @brief Assignment operator from another VecIndex.
         * 
         * @tparam Derived The type of the input vector.
         * @param vec_in The input vector to assign from.
         * @return IndexVecView<Index>& Reference to this object.
         */
        template <typename Derived>
        IndexVecView<Index>& operator=(const VecIndex<Derived, Index> &vec_in)
        {
            fatrop_dbg_assert(this->m() == vec_in.m() &&
                              "Vectors must be same size for assignment");
            for (Index i = 0; i < m(); i++)
            {
                vec_[ai_ + i] = vec_in[i];
            }
            return *this;
        }
        
        /**
         * @brief Assignment operator from another IndexVecView.
         * 
         * @param vec_in The input IndexVecView to assign from.
         * @return IndexVecView<Index>& Reference to this object.
         */
        IndexVecView<Index>& operator=(const IndexVecView& vec_in)
        {
            return *this = static_cast<const VecIndex<IndexVecView<Index>, Index>&>(vec_in);
        }
        
        VecIndexAllocated<Index>& vec_;  ///< Reference to the underlying VecIndexAllocated
        Index ai_;                       ///< Starting index of this view
        Index m_;                        ///< Size of this view
    };

    /**
     * @brief An allocated vector of indices that also provides a view interface.
     * 
     * This class combines the functionality of std::vector with the interface
     * provided by IndexVecView.
     * 
     * @tparam Index The index type used in the vector.
     */
    template <typename Index>
    class VecIndexAllocated : public IndexVecView<Index>,
                              public std::vector<Index>
    {
    public:
        /**
         * @brief Constructs a VecIndexAllocated from an rvalue std::vector.
         * 
         * @param vec The vector to move from.
         */
        VecIndexAllocated(std::vector<Index>&& vec): IndexVecView<Index>(*this, 0, vec.size()), std::vector<Index>(std::move(vec)) {};
        
        /**
         * @brief Constructs a VecIndexAllocated from a const reference to an std::vector.
         * 
         * @param vec The vector to copy from.
         */
        VecIndexAllocated(const std::vector<Index>& vec): IndexVecView<Index>(*this, 0, vec.size()), std::vector<Index>(vec) {};
        
        using std::vector<Index>::operator[];
        
        /**
         * @brief Returns the size of the vector.
         * 
         * @return Index The number of elements in the vector.
         */
        Index m() const { return this->size(); }
    };
    // implementation

    template <typename Index>
        Index IndexVecView<Index>::operator[](const Index i) const { return vec_[i+ai_]; }
    template <typename Index>
        Index IndexVecView<Index>::m() const { return m_; }
} // namespace fatrop

#endif //__fatrop_common_index_vector__
