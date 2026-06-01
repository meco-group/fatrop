#pragma once
#include <vector>
#include <memory>
#include <cassert>
#include <type_traits>

#define CASADI_C_FUNCTIONS(name) \
    casadi_c_wrap::CasadiCFunctions{&name##_n_in, &name##_n_out, &name##_work, \
        &name##_sparsity_in, &name##_sparsity_out, &name##_incref, &name##_decref, \
        &name##_checkout, &name##_release, &name}

namespace casadi_c_wrap
{
    using casadi_int_t = long long int; 
    using casadi_real_t = double;

    // Non-owning view over contiguous data (C++17 compatible)
    template <typename T>
    struct span
    {
        T* data_ = nullptr;
        size_t size_ = 0;

        span() = default;
        span(T* data, size_t size) : data_(data), size_(size) {}
        span(std::vector<std::remove_const_t<T>> &v) : data_(v.data()), size_(v.size()) {}
        span(const std::vector<std::remove_const_t<T>> &v) : data_(const_cast<T*>(v.data())), size_(v.size()) {}

        T* data() const { return data_; }
        size_t size() const { return size_; }
        T& operator[](size_t i) const { return data_[i]; }
    };

    struct CasadiCFunctions
    {
        // n_in
        casadi_int_t (*n_in)();
        // n_out
        casadi_int_t (*n_out)();
        // work
        int (*work)(casadi_int_t *sz_arg, casadi_int_t *sz_res,
                    casadi_int_t *sz_iw, casadi_int_t *sz_w);
        // sparsity_in
        const casadi_int_t *(*sparsity_in)(casadi_int_t i);
        // sparsity_out
        const casadi_int_t *(*sparsity_out)(casadi_int_t i);
        // incref
        void (*incref)();
        // decref
        void (*decref)();
        // checkout
        int (*checkout)();
        // release
        void (*release)(int id);
        // eval
        int (*eval)(const casadi_real_t **arg, casadi_real_t **res, casadi_int_t *iw, casadi_real_t *w, int id);
    };

    class CasadiCFunctionWrap
    {
    public:
        CasadiCFunctionWrap(const CasadiCFunctions &funcs) : funcs_(funcs)
        {
            // Increment reference count for thread-safe memory management
            funcs_.incref();
            mem_id_ = funcs_.checkout();
            allocate_memory();
        }

        ~CasadiCFunctionWrap()
        {
            funcs_.release(mem_id_);
            funcs_.decref();
        }
        void allocate_arg_res(std::vector<std::vector<casadi_real_t>> &arg,
                              std::vector<std::vector<casadi_real_t>> &res) const
        {
            casadi_int_t n_in = funcs_.n_in();
            casadi_int_t n_out = funcs_.n_out();
            arg.resize(n_in);
            for (casadi_int_t i = 0; i < n_in; ++i)
                arg[i].resize(get_nnz(funcs_.sparsity_in(i)), 0.0);
            res.resize(n_out);
            for (casadi_int_t i = 0; i < n_out; ++i)
                res[i].resize(get_nnz(funcs_.sparsity_out(i)));
        }

        template <typename Buffer>
        bool check_arg(const std::vector<Buffer> &arg) const
        {
            casadi_int_t n_in = funcs_.n_in();
            if (static_cast<casadi_int_t>(arg.size()) != n_in) return false;
            for (casadi_int_t i = 0; i < n_in; ++i)
                if (static_cast<casadi_int_t>(arg[i].size()) != get_nnz(funcs_.sparsity_in(i)))
                    return false;
            return true;
        }

        template <typename Buffer>
        bool check_arg_res(const std::vector<Buffer> &arg,
                           const std::vector<Buffer> &res) const
        {
            casadi_int_t n_out = funcs_.n_out();
            if (!check_arg(arg)) return false;
            if (static_cast<casadi_int_t>(res.size()) != n_out) return false;
            for (casadi_int_t i = 0; i < n_out; ++i)
                if (static_cast<casadi_int_t>(res[i].size()) != get_nnz(funcs_.sparsity_out(i)))
                    return false;
            return true;
        }

        template<typename Buffer>
        bool operator()(const std::vector<Buffer> &arg,
                        std::vector<Buffer> &result) 
        {
            assert(check_arg_res(arg, result) && "arg/result sizes do not match function signature");

            for (size_t i = 0; i < arg.size(); ++i)
                memory_.arg[i] = arg[i].data();

            for (size_t i = 0; i < result.size(); ++i)
                memory_.res[i] = result[i].data();

            int eval_result = funcs_.eval(memory_.arg.data(), memory_.res.data(),
                                          memory_.iw.data(), memory_.w.data(), mem_id_);
            return eval_result == 0;
        }

    private:
        static casadi_int_t get_nnz(const casadi_int_t *sp)
        {
            if (!sp)
                return 0;
            casadi_int_t nrow = sp[0];
            casadi_int_t ncol = sp[1];
            return nrow * ncol;
        }

        // Structure to hold all allocated memory for function evaluation
        struct
        {
            std::vector<const casadi_real_t *> arg;
            std::vector<casadi_real_t *> res;
            std::vector<casadi_int_t> iw;
            std::vector<casadi_real_t> w;
        } memory_;
        void allocate_memory()
        {
            casadi_int_t n_in = funcs_.n_in();
            casadi_int_t n_out = funcs_.n_out();

            // Get work vector sizes
            casadi_int_t sz_arg = n_in, sz_res = n_out, sz_iw = 0, sz_w = 0;
            funcs_.work(&sz_arg, &sz_res, &sz_iw, &sz_w);

            // Allocate work vectors
            memory_.arg.resize(sz_arg);
            memory_.res.resize(sz_res);
            memory_.iw.resize(sz_iw);
            memory_.w.resize(sz_w);
        }

        CasadiCFunctions funcs_;
        int mem_id_;
    };

    // Specialization: holds vector buffers for bridging C <-> C++
    template <typename T>
    class EvalData
    {
    public:
        std::vector<std::vector<casadi_real_t>> arg_buffers;
        std::vector<std::vector<casadi_real_t>> res_buffers;
    };

    // note that MemoryManager checkout is not thread-safe
    template <typename T>
    class MemoryManager
    {
    private:
        std::vector<std::unique_ptr<T>> pool_;
        std::vector<bool> in_use_;

        MemoryManager(const MemoryManager &) = delete;
        MemoryManager &operator=(const MemoryManager &) = delete;
        MemoryManager() {
            // pre-create slot 0 for metadata queries
            pool_.push_back(std::make_unique<T>());
            in_use_.push_back(false);
        }

        static MemoryManager& instance_() {
            static MemoryManager manager;
            return manager;
        }

        int checkout_() {
            for (size_t i = 0; i < in_use_.size(); ++i) {
                if (!in_use_[i]) {
                    in_use_[i] = true;
                    return static_cast<int>(i);
                }
            }
            pool_.push_back(std::make_unique<T>());
            in_use_.push_back(true);
            return static_cast<int>(pool_.size() - 1);
        }

        void release_(int id) { in_use_[id] = false; }
        T& get_(int id) { return *pool_[id]; }

    public:
        static int checkout() { return instance_().checkout_(); }
        static void release(int id) { instance_().release_(id); }
        static T& get(int id) { return instance_().get_(id); }
    };

    std::vector<casadi_int_t> dense_sparsity(casadi_int_t nrow, casadi_int_t ncol)
    {
        // CCS format: [nrow, ncol, colind(ncol+1), row(nrow*ncol)]
        std::vector<casadi_int_t> sp(2 + (ncol + 1) + nrow * ncol);
        sp[0] = nrow; sp[1] = ncol;
        for (casadi_int_t i = 0; i <= ncol; ++i) sp[2 + i] = i * nrow;
        for (casadi_int_t c = 0; c < ncol; ++c)
            for (casadi_int_t r = 0; r < nrow; ++r)
                sp[2 + (ncol + 1) + c * nrow + r] = r;
        return sp;
    }

    // macro to expose a C++ class (implementing EvalAbstract) as casadi c functions
    // data_manager: an EvalDataManager<EvalAbstract> that handles checkout/release
    #define EXPOSE_CASADI_C_FUNCTION_WRAP(name, wrap_class) \
        namespace name##_wrap_ns { \
            using memory_type = typename decltype(wrap_class)::memory_type; \
            using data_manager = MemoryManager<memory_type>; \
        } \
        extern "C" { \
        CASADI_SYMBOL_EXPORT long long int name##_n_in() { return wrap_class.n_in(name##_wrap_ns::data_manager::get(0)); } \
        CASADI_SYMBOL_EXPORT long long int name##_n_out() { return wrap_class.n_out(name##_wrap_ns::data_manager::get(0)); } \
        CASADI_SYMBOL_EXPORT int name##_work(long long int *sz_arg, long long int *sz_res, \
            long long int *sz_iw, long long int *sz_w) { \
            *sz_arg = wrap_class.n_in(name##_wrap_ns::data_manager::get(0)); \
            *sz_res = wrap_class.n_out(name##_wrap_ns::data_manager::get(0)); \
            *sz_iw = 0; *sz_w = 0; return 0; } \
        CASADI_SYMBOL_EXPORT const long long int* name##_sparsity_in(long long int i) { return wrap_class.sparsity_in(i, name##_wrap_ns::data_manager::get(0)); } \
        CASADI_SYMBOL_EXPORT const long long int* name##_sparsity_out(long long int i) { return wrap_class.sparsity_out(i, name##_wrap_ns::data_manager::get(0)); } \
        CASADI_SYMBOL_EXPORT void name##_incref() { } \
        CASADI_SYMBOL_EXPORT void name##_decref() { } \
        CASADI_SYMBOL_EXPORT int name##_checkout() { \
            int id = name##_wrap_ns::data_manager::checkout(); \
            wrap_class.init_memory(name##_wrap_ns::data_manager::get(id)); \
            return id; \
        } \
        CASADI_SYMBOL_EXPORT void name##_release(int id) { name##_wrap_ns::data_manager::release(id); } \
        CASADI_SYMBOL_EXPORT int name(const double **arg, double **res, \
            long long int *iw, double *w, int id) { \
            auto &data = name##_wrap_ns::data_manager::get(id); \
            wrap_class(arg, res, data); \
            return 0; \
        } \
    }
} // namespace casadi_c_wrap