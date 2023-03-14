#include <memory>
template<typename T1, typename T2>
void assign_shared_ptr(std::shared_ptr<T1>& lhs, const std::shared_ptr<T2>& rhs){
        lhs = rhs;
}