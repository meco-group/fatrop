#ifndef SMARTPOINTERINCLUDED
#define SMARTPOINTERINCLUDED
namespace fatrop
{
    class RefCountedObj
    {
    public:
        void IncrRef()
        {
            ref_count_++;
        }
        void DecrRef()
        {
            ref_count_--;
        }
        int RefCount()
        {
            return ref_count_;
        }
        virtual ~RefCountedObj(){};

    private:
        int ref_count_ = 0;
    };
    template <typename T>
    class RefCountPtr
    {
    public:
        RefCountPtr()
        {
            ptr_ = NULL;
        }
        RefCountPtr(T *raw_ptr)
        {
            ptr_ = raw_ptr;
            ptr_->IncrRef();
        }
        RefCountPtr(const RefCountPtr<T> &copy) : ptr_(copy)
        {
            ptr_->IncrRef();
        }
        inline T *operator->() const
        {
            return ptr_;
        }
        inline operator T *() const
        {
            return ptr_;
        }
        inline T &operator*() const
        {
            return *ptr_;
        }
        T *GetRawPtr() const
        {
            return ptr_;
        }
        RefCountPtr<T> &operator=(const RefCountPtr<T> &other)
        {
            RemoveRef();
            this->ptr_ = other.GetRawPtr();
            this->ptr_->IncrRef();
            return *this;
        }
        void RemoveRef()
        {
            if (this->ptr_ != NULL)
            {
                ptr_->DecrRef();
                if (ptr_->RefCount() == 0)
                {
                    delete ptr_;
                }
            }
        }
        ~RefCountPtr()
        {
            RemoveRef();
        }

    private:
        T *ptr_ = nullptr;
    };
} // namespace fatrop
#endif // SMARTPOINTERINCLUDED