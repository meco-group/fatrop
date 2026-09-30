#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <vector>

namespace fatrop_benchmark::ad
{
constexpr int capacity = 16;

// Each OCP uses only as many forward directions as its stage variables.
template <int D> struct Forward
{
    double value = 0.0;
    std::array<double, D> tangent{};

    Forward() = default;
    Forward(double value) : value(value) {}
    static Forward variable(double value, int index)
    {
        Forward result(value);
        result.tangent[index] = 1.0;
        return result;
    }
};

template <int D> Forward<D> operator+(const Forward<D> &a, const Forward<D> &b)
{
    Forward<D> r(a.value + b.value);
    for (int i = 0; i < D; ++i) r.tangent[i] = a.tangent[i] + b.tangent[i];
    return r;
}
template <int D> Forward<D> operator-(const Forward<D> &a)
{
    Forward<D> r(-a.value);
    for (int i = 0; i < D; ++i) r.tangent[i] = -a.tangent[i];
    return r;
}
template <int D> Forward<D> operator-(const Forward<D> &a, const Forward<D> &b) { return a + (-b); }
template <int D> Forward<D> operator*(const Forward<D> &a, const Forward<D> &b)
{
    Forward<D> r(a.value * b.value);
    for (int i = 0; i < D; ++i)
        r.tangent[i] = a.tangent[i] * b.value + a.value * b.tangent[i];
    return r;
}
template <int D> Forward<D> operator/(const Forward<D> &a, const Forward<D> &b)
{
    Forward<D> r(a.value / b.value);
    for (int i = 0; i < D; ++i)
        r.tangent[i] = (a.tangent[i] - r.value * b.tangent[i]) / b.value;
    return r;
}
template <int D> Forward<D> operator+(const Forward<D> &a, double b) { return a + Forward<D>(b); }
template <int D> Forward<D> operator+(double a, const Forward<D> &b) { return Forward<D>(a) + b; }
template <int D> Forward<D> operator-(const Forward<D> &a, double b) { return a - Forward<D>(b); }
template <int D> Forward<D> operator-(double a, const Forward<D> &b) { return Forward<D>(a) - b; }
template <int D> Forward<D> operator*(const Forward<D> &a, double b) { return a * Forward<D>(b); }
template <int D> Forward<D> operator*(double a, const Forward<D> &b) { return Forward<D>(a) * b; }
template <int D> Forward<D> operator/(const Forward<D> &a, double b) { return a / Forward<D>(b); }
template <int D> Forward<D> operator/(double a, const Forward<D> &b) { return Forward<D>(a) / b; }

template <class Scalar> class Tape;

template <class Scalar> struct Number
{
    Scalar value{};
    Tape<Scalar> *tape = nullptr;
    int node = -1;

    Number() = default;
    Number(double value) : value(value) {}
    Number(Scalar value, Tape<Scalar> *tape, int node)
        : value(value), tape(tape), node(node) {}
};

template <class Scalar> class Tape
{
    struct Node
    {
        int left = -1, right = -1;
        Scalar left_weight{}, right_weight{}, adjoint{};
    };
    std::vector<Node> nodes_;

public:
    Tape() { nodes_.reserve(1024); }
    void reset() { nodes_.clear(); }

    Number<Scalar> variable(Scalar value)
    {
        const int index = int(nodes_.size());
        nodes_.emplace_back();
        return {value, this, index};
    }

    Number<Scalar> record(Scalar value, const Number<Scalar> &a,
                          const Number<Scalar> &b, Scalar da, Scalar db)
    {
        assert(!a.tape || !b.tape || a.tape == b.tape);
        Tape *owner = a.tape ? a.tape : b.tape;
        assert(owner == this);
        const int index = int(owner->nodes_.size());
        owner->nodes_.push_back({a.node, b.node, da, db, Scalar(0.0)});
        return {value, owner, index};
    }

    void backward(const Number<Scalar> &output)
    {
        for (auto &node : nodes_) node.adjoint = Scalar(0.0);
        if (!output.tape) return;
        assert(output.tape == this);
        nodes_[output.node].adjoint = Scalar(1.0);
        for (int i = int(nodes_.size()) - 1; i >= 0; --i)
        {
            const auto &node = nodes_[i];
            if (node.left >= 0)
                nodes_[node.left].adjoint = nodes_[node.left].adjoint
                    + node.adjoint * node.left_weight;
            if (node.right >= 0)
                nodes_[node.right].adjoint = nodes_[node.right].adjoint
                    + node.adjoint * node.right_weight;
        }
    }

    Scalar adjoint(const Number<Scalar> &variable) const
    { return nodes_[variable.node].adjoint; }
};

template <class S> Number<S> record(S value, const Number<S> &a,
                                    const Number<S> &b, S da, S db)
{
    Tape<S> *tape = a.tape ? a.tape : b.tape;
    return tape ? tape->record(value, a, b, da, db) : Number<S>{value, nullptr, -1};
}
template <class S> Number<S> operator+(const Number<S> &a, const Number<S> &b)
{ return record(a.value + b.value, a, b, S(1.0), S(1.0)); }
template <class S> Number<S> operator-(const Number<S> &a, const Number<S> &b)
{ return record(a.value - b.value, a, b, S(1.0), S(-1.0)); }
template <class S> Number<S> operator-(const Number<S> &a)
{ return Number<S>(0.0) - a; }
template <class S> Number<S> operator*(const Number<S> &a, const Number<S> &b)
{ return record(a.value * b.value, a, b, b.value, a.value); }
template <class S> Number<S> operator/(const Number<S> &a, const Number<S> &b)
{
    const S inverse = S(1.0) / b.value;
    return record(a.value * inverse, a, b, inverse, -a.value * inverse * inverse);
}
template <class S> Number<S> operator+(double a, const Number<S> &b) { return Number<S>(a) + b; }
template <class S> Number<S> operator+(const Number<S> &a, double b) { return a + Number<S>(b); }
template <class S> Number<S> operator-(double a, const Number<S> &b) { return Number<S>(a) - b; }
template <class S> Number<S> operator-(const Number<S> &a, double b) { return a - Number<S>(b); }
template <class S> Number<S> operator*(double a, const Number<S> &b) { return Number<S>(a) * b; }
template <class S> Number<S> operator*(const Number<S> &a, double b) { return a * Number<S>(b); }
template <class S> Number<S> operator/(double a, const Number<S> &b) { return Number<S>(a) / b; }
template <class S> Number<S> operator/(const Number<S> &a, double b) { return a / Number<S>(b); }

template <class S> Number<S> unary(const Number<S> &a, S value, S derivative)
{ return record(value, a, Number<S>{}, derivative, S(0.0)); }
template <class S> Number<S> sin(const Number<S> &a)
{
    using std::sin; using std::cos;
    return unary(a, sin(a.value), cos(a.value));
}
template <class S> Number<S> cos(const Number<S> &a)
{
    using std::sin; using std::cos;
    return unary(a, cos(a.value), -sin(a.value));
}
template <class S> Number<S> exp(const Number<S> &a)
{
    using std::exp;
    const S value = exp(a.value);
    return unary(a, value, value);
}
template <class S> Number<S> tan(const Number<S> &a)
{
    using std::tan;
    const S value = tan(a.value);
    return unary(a, value, S(1.0) + value * value);
}
template <class S> Number<S> tanh(const Number<S> &a)
{
    using std::tanh;
    const S value = tanh(a.value);
    return unary(a, value, S(1.0) - value * value);
}

template <int D> Forward<D> forward_unary(const Forward<D> &a, double value, double derivative)
{
    Forward<D> r(value);
    for (int i = 0; i < D; ++i) r.tangent[i] = derivative * a.tangent[i];
    return r;
}
template <int D> Forward<D> sin(const Forward<D> &a)
{ return forward_unary(a, std::sin(a.value), std::cos(a.value)); }
template <int D> Forward<D> cos(const Forward<D> &a)
{ return forward_unary(a, std::cos(a.value), -std::sin(a.value)); }
template <int D> Forward<D> exp(const Forward<D> &a)
{ const double e = std::exp(a.value); return forward_unary(a, e, e); }
template <int D> Forward<D> tan(const Forward<D> &a)
{
    const double t = std::tan(a.value);
    return forward_unary(a, t, 1.0 + t * t);
}
template <int D> Forward<D> tanh(const Forward<D> &a)
{
    const double t = std::tanh(a.value);
    return forward_unary(a, t, 1.0 - t * t);
}
inline double sin(double a) { return std::sin(a); }
inline double cos(double a) { return std::cos(a); }
inline double exp(double a) { return std::exp(a); }
inline double tan(double a) { return std::tan(a); }
inline double tanh(double a) { return std::tanh(a); }

class HessianTape;

struct HessianNumber
{
    double value = 0.0;
    HessianTape *tape = nullptr;
    int node = -1;

    HessianNumber() = default;
    HessianNumber(double value) : value(value) {}
    HessianNumber(double value, HessianTape *tape, int node)
        : value(value), tape(tape), node(node) {}
};

// Forward-over-adjoint Hessian, one direction at a time.
class HessianTape
{
    struct Node
    {
        int left = -1, right = -1;
        double dl = 0.0, dr = 0.0;
        double dll = 0.0, dlr = 0.0, drr = 0.0;
        double adjoint = 0.0;
    };
    std::vector<Node> nodes_;
    std::vector<double> tangent_, adjoint_tangent_;

public:
    HessianTape() { nodes_.reserve(256); }
    void reset() { nodes_.clear(); }
    HessianNumber variable(double value)
    {
        const int index = int(nodes_.size());
        nodes_.emplace_back();
        return {value, this, index};
    }
    HessianNumber record(double value, const HessianNumber &a,
                         const HessianNumber &b, double dl, double dr,
                         double dll, double dlr, double drr)
    {
        assert(!a.tape || !b.tape || a.tape == b.tape);
        assert((a.tape ? a.tape : b.tape) == this);
        const int index = int(nodes_.size());
        nodes_.push_back({a.node, b.node, dl, dr, dll, dlr, drr, 0.0});
        return {value, this, index};
    }
    std::array<double, capacity * capacity> hessian(
        const HessianNumber &output, const int *variables, int dimension)
    {
        std::array<double, capacity * capacity> result{};
        if (!output.tape) return result;
        assert(output.tape == this);
        for (auto &node : nodes_) node.adjoint = 0.0;
        nodes_[output.node].adjoint = 1.0;
        for (int i = int(nodes_.size()) - 1; i >= 0; --i)
        {
            const Node &node = nodes_[i];
            if (node.left >= 0) nodes_[node.left].adjoint += node.adjoint * node.dl;
            if (node.right >= 0) nodes_[node.right].adjoint += node.adjoint * node.dr;
        }
        tangent_.resize(nodes_.size());
        adjoint_tangent_.resize(nodes_.size());
        for (int column = 0; column < dimension; ++column)
        {
            std::fill(tangent_.begin(), tangent_.end(), 0.0);
            tangent_[variables[column]] = 1.0;
            for (int i = 0; i < int(nodes_.size()); ++i)
            {
                const Node &node = nodes_[i];
                if (node.left < 0 && node.right < 0) continue;
                tangent_[i] = (node.left >= 0 ? node.dl * tangent_[node.left] : 0.0)
                            + (node.right >= 0 ? node.dr * tangent_[node.right] : 0.0);
            }
            std::fill(adjoint_tangent_.begin(), adjoint_tangent_.end(), 0.0);
            for (int i = int(nodes_.size()) - 1; i >= 0; --i)
            {
                const Node &node = nodes_[i];
                const double left_tangent = node.left >= 0 ? tangent_[node.left] : 0.0;
                const double right_tangent = node.right >= 0 ? tangent_[node.right] : 0.0;
                if (node.left >= 0)
                    adjoint_tangent_[node.left] += adjoint_tangent_[i] * node.dl
                        + node.adjoint * (node.dll * left_tangent + node.dlr * right_tangent);
                if (node.right >= 0)
                    adjoint_tangent_[node.right] += adjoint_tangent_[i] * node.dr
                        + node.adjoint * (node.dlr * left_tangent + node.drr * right_tangent);
            }
            for (int row = 0; row < dimension; ++row)
                result[row * capacity + column] = adjoint_tangent_[variables[row]];
        }
        return result;
    }
};

inline HessianNumber hessian_record(double value, const HessianNumber &a,
                                    const HessianNumber &b, double dl, double dr,
                                    double dll, double dlr, double drr)
{
    HessianTape *tape = a.tape ? a.tape : b.tape;
    return tape ? tape->record(value, a, b, dl, dr, dll, dlr, drr)
                : HessianNumber(value);
}
inline HessianNumber operator+(const HessianNumber &a, const HessianNumber &b)
{ return hessian_record(a.value + b.value, a, b, 1.0, 1.0, 0.0, 0.0, 0.0); }
inline HessianNumber operator-(const HessianNumber &a, const HessianNumber &b)
{ return hessian_record(a.value - b.value, a, b, 1.0, -1.0, 0.0, 0.0, 0.0); }
inline HessianNumber operator-(const HessianNumber &a)
{ return HessianNumber(0.0) - a; }
inline HessianNumber operator*(const HessianNumber &a, const HessianNumber &b)
{ return hessian_record(a.value * b.value, a, b, b.value, a.value, 0.0, 1.0, 0.0); }
inline HessianNumber operator/(const HessianNumber &a, const HessianNumber &b)
{
    const double inv = 1.0 / b.value;
    return hessian_record(a.value * inv, a, b, inv, -a.value * inv * inv,
                          0.0, -inv * inv, 2.0 * a.value * inv * inv * inv);
}
inline HessianNumber operator+(double a, const HessianNumber &b) { return HessianNumber(a) + b; }
inline HessianNumber operator+(const HessianNumber &a, double b) { return a + HessianNumber(b); }
inline HessianNumber operator-(double a, const HessianNumber &b) { return HessianNumber(a) - b; }
inline HessianNumber operator-(const HessianNumber &a, double b) { return a - HessianNumber(b); }
inline HessianNumber operator*(double a, const HessianNumber &b) { return HessianNumber(a) * b; }
inline HessianNumber operator*(const HessianNumber &a, double b) { return a * HessianNumber(b); }
inline HessianNumber operator/(double a, const HessianNumber &b) { return HessianNumber(a) / b; }
inline HessianNumber operator/(const HessianNumber &a, double b) { return a / HessianNumber(b); }

inline HessianNumber hessian_unary(const HessianNumber &a, double value,
                                   double first, double second)
{ return hessian_record(value, a, HessianNumber{}, first, 0.0, second, 0.0, 0.0); }
inline HessianNumber sin(const HessianNumber &a)
{ return hessian_unary(a, std::sin(a.value), std::cos(a.value), -std::sin(a.value)); }
inline HessianNumber cos(const HessianNumber &a)
{ return hessian_unary(a, std::cos(a.value), -std::sin(a.value), -std::cos(a.value)); }
inline HessianNumber exp(const HessianNumber &a)
{ const double e = std::exp(a.value); return hessian_unary(a, e, e, e); }
inline HessianNumber tan(const HessianNumber &a)
{
    const double t = std::tan(a.value), first = 1.0 + t * t;
    return hessian_unary(a, t, first, 2.0 * t * first);
}
inline HessianNumber tanh(const HessianNumber &a)
{
    const double t = std::tanh(a.value), first = 1.0 - t * t;
    return hessian_unary(a, t, first, -2.0 * t * first);
}
}
