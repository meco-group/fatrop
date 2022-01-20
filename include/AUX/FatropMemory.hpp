/**
 *  @file FatropMemory.hpp
 *  Memory allocation of fatrop. 
 *
*/
#ifndef FATROP_MEMORY_INCLUDED
#define FATROP_MEMORY_INCLUDED
#include <iostream>
#include <vector>
#include "FatropVector.hpp"
#include <utility>
using namespace std;
namespace fatrop
{
    // /** \brief fatrop memory element, used for allocation of memory*/
    // class MemoryElBase
    // {
    // public:
    //     /** \brief calculate memory size*/
    //     virtual int memory_size() const = 0;
    //     /** \brief set up memory element and advance pointer */
    //     virtual void set_up(char *&data_p) = 0;
    // };

    // /** \brief allocates and initializes fatrop_memory_el's in one big chunk of memory */
    // class MemoryAllocator
    // {
    // public:
    //     /// add a fatrop_memory_el to allocator
    //     void add(MemoryElBase &fme) { fatrop_memory_el_vector.push_back(&fme); };
    //     /// calculate memory size needed to allocate all fatrop_memroy_el's of this allocator
    //     unsigned long int memory_size()
    //     {
    //         int size = fatrop_memory_el_vector.size();
    //         unsigned long int res = 0;
    //         for (int i = 0; i < size; i++)
    //         {
    //             res += (unsigned long int) fatrop_memory_el_vector.at(i)->memory_size();
    //         }
    //         return res;
    //     }
    //     /// allocate memory and set up all fatrop_memory_el's of this container
    //     void allocate()
    //     {
    //         int size = fatrop_memory_el_vector.size();
    //         unsigned long int mem_size = this->memory_size();
    //         cout << "allocating buffer of size " << mem_size / 1024 << " kibytes and " << mem_size % 1024 << " byte" << endl;
    //         // storage_mem = malloc(mem_size);
    //         // char *c_ptr = (char *)storage_mem;
    //         for (int i = 0; i < size; i++)
    //         {
    //             storage_mem = malloc(fatrop_memory_el_vector.at(i) ->memory_size());
    //             char *c_ptr = (char *)storage_mem;
    //             fatrop_memory_el_vector.at(i)->set_up(c_ptr);
    //         }
    //     };

    // private:
    //     /// stores which fatrop_memory_el's are contained
    //     vector<MemoryElBase *> fatrop_memory_el_vector;
    //     /// pointer to memory location of allocated data
    //     void *storage_mem = NULL;

    // public:
    //     /// destructor
    //     ~MemoryAllocator()
    //     {
    //         free(storage_mem);
    //     }
    // };

    // /** \brief fatrop memory element template */
    // /// example usage:
    // /// fatrop_memory_el<int> test(5);
    // /// represents five integers in memory
    // ///
    // /// default constructor should be available
    // template <typename T>
    // class FatropMemoryEl : public MemoryElBase
    // {
    // public:
    //     // E&& is a universal reference
    //     template<typename E>
    //     FatropMemoryEl<T>(int size, E&&init_values, MemoryAllocator &fma) : size(size), init_values_(forward<E>(init_values))
    //     {
    //         fma.add(*this);
    //     }
    //     int memory_size() const
    //     {
    //         return size * sizeof(T);
    //     }
    //     void set_up(char *&data_p)
    //     {
    //         data = (T *)data_p;
    //         for (int i = 0; i < size; i++)
    //         {
    //             T *t_p = (T *)data_p;
    //             *t_p = init_values_.at(i); // default constructor
    //             data_p += sizeof(T);       // advance pointer
    //         }
    //     }
    //     explicit operator T *() const { return data; };

    // private:
    //     const int size;
    //     FatropVector<T> init_values_;
    //     T *data = NULL;
    // };
} // namespace fatrop
#endif //FATROP_MEMORY_INCLUDED