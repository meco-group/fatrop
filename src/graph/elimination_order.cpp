//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/elimination_order.hpp"

#include "fatrop/graph/block_sparsity.hpp"

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

using namespace fatrop;

// Approximate minimum degree (Amestoy, Davis, Duff, SIAM J. Matrix Anal. Appl.
// 17(4), 1996) on the block graph, with every degree weighted by block size.
//
// The elimination graph is kept as a quotient graph: an eliminated (super)
// variable p becomes an element whose variable list L_p stands for the clique
// it created, so fill is never formed explicitly. A variable i keeps a list of
// adjacent elements E_i and of adjacent variables A_i; its external degree is
// bounded by
//   d_i <= min(rem - w_i, d_i_old + |L_p \ i|,
//              |A_i| + |L_p \ i| + sum_{e in E_i \ p} |L_e \ L_p|),
// where |.| sums the block sizes of the (super)variables in the set and rem is
// the remaining weight. Elements covered by L_p are absorbed, variables only
// adjacent to p are eliminated together with p (mass elimination) and
// variables with identical adjacency are merged into supervariables. Degree
// buckets give O(1) pivot selection; the bucket heads are LIFO, which makes
// the result deterministic.
//
// The members of a supervariable are emitted consecutively, so they form a
// chain in the elimination tree with identical column patterns (a supernode).
namespace
{
enum NodeState : unsigned char
{
    VARIABLE, // principal (super)variable
    ELEMENT,  // eliminated, not yet absorbed
    DEAD      // absorbed element or non-principal variable
};
} // namespace

BlockEliminationOrder::BlockEliminationOrder(const BlockSparsityPattern &sp)
    : perm_(sp.num_blocks()), inv_perm_(sp.num_blocks())
{
    const Index n = sp.num_blocks();
    if (n == 0)
        return;
    // var[i] / elem[i]: A_i / E_i of a variable; var[e] holds L_e of an element.
    std::vector<std::vector<Index>> var(n), elem(n);
    std::vector<NodeState> state(n, VARIABLE);
    // w: scalar size, nv: number of blocks of a supervariable (for elements:
    // summed over L_e; both are invariant until the element is absorbed).
    std::vector<Index> w(sp.block_sizes()), nv(n, 1), deg(n, 0);
    std::vector<Index> el_w(n, 0), el_nv(n, 0);
    // Supervariable members as linked lists.
    std::vector<Index> member_next(n, -1), member_last(n);
    // Degree buckets (doubly linked lists).
    const Index max_deg = sp.total_size();
    std::vector<Index> head(max_deg + 1, -1), bnext(n, -1), bprev(n, -1);
    // Marks and |L_e \ L_p| (weighted and in blocks).
    std::vector<Index> mark(n, -1), ext_mark(n, -1), ext_w(n, 0), ext_nv(n, 0);
    std::vector<Index> deg_new(n, 0);
    std::vector<Index> lp;
    std::vector<std::pair<std::uint64_t, Index>> hashed;
    lp.reserve(n);
    hashed.reserve(n);

    Index mindeg = max_deg;
    auto bucket_insert = [&](Index i)
    {
        const Index d = deg[i];
        bnext[i] = head[d];
        bprev[i] = -1;
        if (head[d] != -1)
            bprev[head[d]] = i;
        head[d] = i;
        mindeg = std::min(mindeg, d);
    };
    auto bucket_remove = [&](Index i)
    {
        if (bprev[i] != -1)
            bnext[bprev[i]] = bnext[i];
        else
            head[deg[i]] = bnext[i];
        if (bnext[i] != -1)
            bprev[bnext[i]] = bprev[i];
    };
    Index n_out = 0;
    auto emit = [&](Index i)
    {
        for (Index k = i; k != -1; k = member_next[k])
        {
            perm_[n_out] = k;
            inv_perm_[k] = n_out++;
        }
    };
    auto release = [](std::vector<Index> &v) { std::vector<Index>().swap(v); };

    for (Index i = 0; i < n; ++i)
    {
        var[i] = sp.neighbors(i);
        member_last[i] = i;
        for (Index j : var[i])
            deg[i] += w[j];
    }
    // Insert in reverse so that, among equal degrees, low indices come first.
    for (Index i = n - 1; i >= 0; --i)
        bucket_insert(i);

    Index rem = sp.total_size();
    Index stamp = 0;
    while (n_out < n)
    {
        while (head[mindeg] == -1)
            ++mindeg;
        const Index p = head[mindeg];
        bucket_remove(p);

        // L_p = (A_p u union_{e in E_p} L_e) \ p; the elements in E_p are absorbed.
        ++stamp;
        lp.clear();
        mark[p] = stamp;
        auto add_to_lp = [&](Index i)
        {
            if (state[i] == VARIABLE && mark[i] != stamp)
            {
                mark[i] = stamp;
                lp.push_back(i);
            }
        };
        for (Index e : elem[p])
        {
            if (state[e] != ELEMENT)
                continue;
            for (Index i : var[e])
                add_to_lp(i);
            state[e] = DEAD;
            release(var[e]);
        }
        for (Index i : var[p])
            add_to_lp(i);
        release(var[p]);
        release(elem[p]);
        state[p] = ELEMENT;
        rem -= w[p];
        emit(p);
        for (Index i : lp)
            bucket_remove(i);

        // |L_e \ L_p| for every element adjacent to L_p.
        for (Index i : lp)
            for (Index e : elem[i])
            {
                if (state[e] != ELEMENT)
                    continue;
                if (ext_mark[e] != stamp)
                {
                    ext_mark[e] = stamp;
                    ext_w[e] = el_w[e];
                    ext_nv[e] = el_nv[e];
                }
                ext_w[e] -= w[i];
                ext_nv[e] -= nv[i];
            }

        // Prune E_i and A_i, absorb elements covered by L_p and eliminate the
        // variables that are adjacent to p only.
        Index lp_w = 0;
        Index lp_nv = 0;
        size_t n_lp = 0;
        for (Index i : lp)
        {
            Index d = 0;
            size_t m = 0;
            auto &E = elem[i];
            for (Index e : E)
            {
                if (state[e] != ELEMENT)
                    continue;
                if (ext_nv[e] == 0)
                {
                    state[e] = DEAD;
                    release(var[e]);
                    continue;
                }
                E[m++] = e;
                d += ext_w[e];
            }
            E.resize(m);
            m = 0;
            auto &A = var[i];
            for (Index j : A)
            {
                if (state[j] != VARIABLE || mark[j] == stamp)
                    continue;
                A[m++] = j;
                d += w[j];
            }
            A.resize(m);
            if (E.empty() && A.empty())
            {
                state[i] = DEAD;
                rem -= w[i];
                release(E);
                release(A);
                emit(i);
                continue;
            }
            deg_new[i] = d;
            lp_w += w[i];
            lp_nv += nv[i];
            lp[n_lp++] = i;
        }
        lp.resize(n_lp);

        for (Index i : lp)
        {
            const Index ext = lp_w - w[i];
            deg[i] = std::min({rem - w[i], deg[i] + ext, deg_new[i] + ext});
            elem[i].push_back(p);
        }

        // Supervariable detection: hash the adjacency, compare within equal hashes.
        hashed.clear();
        for (Index i : lp)
        {
            std::uint64_t h = 0;
            for (Index e : elem[i])
                h += static_cast<std::uint64_t>(e);
            for (Index j : var[i])
                h += static_cast<std::uint64_t>(j);
            hashed.emplace_back(h, i);
        }
        std::sort(hashed.begin(), hashed.end());
        for (size_t a = 0; a < hashed.size(); ++a)
        {
            const Index i = hashed[a].second;
            if (state[i] != VARIABLE)
                continue;
            bool marked = false;
            for (size_t b = a + 1; b < hashed.size() && hashed[b].first == hashed[a].first; ++b)
            {
                const Index j = hashed[b].second;
                if (state[j] != VARIABLE || elem[j].size() != elem[i].size() ||
                    var[j].size() != var[i].size())
                    continue;
                if (!marked)
                {
                    ++stamp;
                    for (Index e : elem[i])
                        mark[e] = stamp;
                    for (Index k : var[i])
                        mark[k] = stamp;
                    marked = true;
                }
                bool same = true;
                for (Index e : elem[j])
                    same = same && mark[e] == stamp;
                for (Index k : var[j])
                    same = same && mark[k] == stamp;
                if (!same)
                    continue;
                // Merge j into i.
                state[j] = DEAD;
                w[i] += w[j];
                nv[i] += nv[j];
                deg[i] = std::max<Index>(deg[i] - w[j], 0);
                member_next[member_last[i]] = j;
                member_last[i] = member_last[j];
                release(elem[j]);
                release(var[j]);
            }
        }

        // L_p keeps the principal variables only.
        n_lp = 0;
        for (Index i : lp)
            if (state[i] == VARIABLE)
            {
                lp[n_lp++] = i;
                bucket_insert(i);
            }
        lp.resize(n_lp);
        var[p].assign(lp.begin(), lp.end());
        el_w[p] = lp_w;
        el_nv[p] = lp_nv;
    }
}
