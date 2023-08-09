#pragma once
namespace fatrop
{
    namespace fatrop_casadi
    {
        template <typename T, typename Derived>
        struct SharedObj: public std::shared_ptr<T>
        {
            SharedObj(const std::shared_ptr<T> &other) : std::shared_ptr<T>(other) {};
            SharedObj() : std::shared_ptr<T>(nullptr) {};
            template <class... Args>
            static Derived create (const Args &...args) { return Derived(std::make_shared<T>(args...)); }
        };
    }; // namespace casadi
}; // namespace fatrop