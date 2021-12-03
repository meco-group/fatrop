/**
 *  @file FatropMemory.hpp
 *  Memory allocation of fatrop. 
 *
*/
#ifndef FATROP_MEMORY_INCLUDED
#define FATROP_MEMORY_INCLUDED
#include <iostream>
#include <vector>
#include <FatropVector.hpp>
using namespace std;
namespace fatrop
{
    /** \brief fatrop memory element, used for allocation of memory*/
    class fatrop_memory_el_base
    {
    public:
        /** \brief calculate memory size*/
        virtual int memory_size() const = 0;
        /** \brief set up memory element and advance pointer */
        virtual void set_up(char *&data_p) = 0;
    };

    /** \brief allocates and initializes fatrop_memory_el's in one big chunk of memory */
    class fatrop_memory_allocator
    {
    public:
        /// add a fatrop_memory_el to allocator
        void add(fatrop_memory_el_base &fme) { fatrop_memory_el_vector.push_back(&fme); };
        /// calculate memory size needed to allocate all fatrop_memroy_el's of this allocator
        int memory_size()
        {
            int size = fatrop_memory_el_vector.size();
            int res = 0;
            for (int i = 0; i < size; i++)
            {
                res += fatrop_memory_el_vector.at(i)->memory_size();
            }
            return res;
        }
        /// allocate memory and set up all fatrop_memory_el's of this container
        void allocate()
        {
            int size = fatrop_memory_el_vector.size();
            int mem_size = this->memory_size();
            cout << "allocating buffer of size " << mem_size / 1024 << " kibytes and " << mem_size % 1024 << " byte" << endl;
            storage_mem = malloc(mem_size);
            char *c_ptr = (char *)storage_mem;
            for (int i = 0; i < size; i++)
            {
                fatrop_memory_el_vector.at(i)->set_up(c_ptr);
            }
        };

    private:
        /// stores which fatrop_memory_el's are contained
        vector<fatrop_memory_el_base *> fatrop_memory_el_vector;
        /// pointer to memory location of allocated data
        void *storage_mem = NULL;

    public:
        /// destructor
        ~fatrop_memory_allocator()
        {
            free(storage_mem);
        }
    };

    /** \brief fatrop memory element template */
    /// example usage:
    /// fatrop_memory_el<int> test(5);
    /// represents five integers in memory
    ///
    /// default constructor should be available
    template <typename T>
    class fatrop_memory_el : public fatrop_memory_el_base
    {
    public:
        template<typename E>
        fatrop_memory_el<T>(int size, const VecExpr<E,T> &init_values, fatrop_memory_allocator &fma) : size(size), init_values_(init_values)
        {
            fma.add(*this);
        }
        fatrop_memory_el<T>(int size, const vector<T> &&init_values, fatrop_memory_allocator &fma) : size(size), init_values_(move(init_values))
        {
            fma.add(*this);
        }
        fatrop_memory_el<T>(int size, const vector<T> &init_values, fatrop_memory_allocator &fma) : size(size), init_values_(init_values)
        {
            fma.add(*this);
        }
        int memory_size() const
        {
            return size * sizeof(T);
        }
        void set_up(char *&data_p)
        {
            data = (T *)data_p;
            for (int i = 0; i < size; i++)
            {
                T *t_p = (T *)data_p;
                *t_p = init_values_.at(i); // default constructor
                data_p += sizeof(T);       // advance pointer
            }
        }
        explicit operator T *() const { return data; };

    private:
        const int size;
        FatropVector<T> init_values_;
        T *data = NULL;
    };
} // namespace fatrop
#endif //FATROP_MEMORY_INCLUDED