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
// buckets (LIFO) give O(1) pivot selection.
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

// sum_{c = a + 1}^{b} c^2: flops of scalar columns with counts a + 1, ..., b.
double sum_squares(double a, double b)
{
    auto s = [](double x) { return x * (x + 1.) * (2. * x + 1.) / 6.; };
    return s(b) - s(a);
}

// One AMD run. Ties between equal degrees are decided by the bucket order,
// which the two flags flip: the initial insertion order and the order in which
// L_p is reinserted. Returns the flops (sum of squared scalar column counts of
// L), which are exact since L_p is the exact pattern of the pivot's columns,
// and the scalar off-diagonal count of every block column of L in count.
double amd(const BlockSparsityPattern &sp, bool ascending_init, bool reverse_reinsert,
           std::vector<Index> &perm, std::vector<Index> &inv_perm, std::vector<Index> &count)
{
    const Index n = sp.num_blocks();
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
            perm[n_out] = k;
            inv_perm[k] = n_out++;
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
    for (Index k = 0; k < n; ++k)
        bucket_insert(ascending_init ? k : n - 1 - k);

    Index rem = sp.total_size();
    Index stamp = 0;
    double flops = 0.;
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
        const Index rem_before = rem;
        rem -= w[p];
        const Index n_first = n_out;
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
                lp[n_lp++] = i;
        lp.resize(n_lp);
        for (size_t q = 0; q < n_lp; ++q)
            bucket_insert(lp[reverse_reinsert ? n_lp - 1 - q : q]);
        var[p].assign(lp.begin(), lp.end());
        el_w[p] = lp_w;
        el_nv[p] = lp_nv;
        // p and the mass-eliminated variables form one chain above L_p.
        flops += sum_squares(lp_w, lp_w + rem_before - rem);
        for (Index q = n_out - 1, below = lp_w; q >= n_first; --q)
        {
            count[perm[q]] = below;
            below += sp.block_size(perm[q]);
        }
    }
    return flops;
}
} // namespace

BlockEliminationOrder::BlockEliminationOrder(const BlockSparsityPattern &sp)
    : perm_(sp.num_blocks()), inv_perm_(sp.num_blocks())
{
    const Index n = sp.num_blocks();
    // AMD is sensitive to tie-breaking on regular meshes; keep the cheapest of
    // four tie-breaking variants.
    std::vector<Index> perm(n), inv_perm(n), count(n), best_count(n);
    double best = -1.;
    for (int variant = 0; variant < 4; ++variant)
    {
        const double flops = amd(sp, variant & 1, variant & 2, perm, inv_perm, count);
        if (best < 0. || flops < best)
        {
            best = flops;
            perm_.swap(perm);
            inv_perm_.swap(inv_perm);
            best_count.swap(count);
        }
    }

    // Postorder the elimination tree (same fill): subtrees become contiguous,
    // so a column directly follows its last child and chains can be merged
    // into supernodes. The child with the largest column count goes last; its
    // pattern overlaps most with the parent's.
    // Elimination tree (Liu) in AMD positions.
    std::vector<Index> parent(n, -1), ancestor(n, -1);
    for (Index k = 0; k < n; ++k)
        for (Index u : sp.neighbors(perm_[k]))
        {
            Index r = inv_perm_[u];
            if (r > k)
                continue;
            while (ancestor[r] != -1 && ancestor[r] != k)
            {
                const Index t = ancestor[r];
                ancestor[r] = k;
                r = t;
            }
            if (ancestor[r] == -1)
            {
                ancestor[r] = k;
                parent[r] = k;
            }
        }
    // Children sorted by increasing column count, as CSR.
    std::vector<Index> child_start(n + 1, 0), child(n), by_count(n);
    for (Index k = 0; k < n; ++k)
        if (parent[k] != -1)
            ++child_start[parent[k] + 1];
    for (Index k = 0; k < n; ++k)
        child_start[k + 1] += child_start[k];
    for (Index k = 0; k < n; ++k)
        by_count[k] = k;
    std::stable_sort(by_count.begin(), by_count.end(), [&](Index a, Index b)
                     { return best_count[perm_[a]] < best_count[perm_[b]]; });
    std::vector<Index> fill(child_start.begin(), child_start.end() - 1);
    for (Index k : by_count)
        if (parent[k] != -1)
            child[fill[parent[k]]++] = k;
    // Iterative depth-first search from the roots.
    std::vector<Index> post, stack, next(n);
    post.reserve(n);
    stack.reserve(n);
    for (Index root = 0; root < n; ++root)
    {
        if (parent[root] != -1)
            continue;
        stack.push_back(root);
        next[root] = child_start[root];
        while (!stack.empty())
        {
            const Index k = stack.back();
            if (next[k] < child_start[k + 1])
            {
                const Index c = child[next[k]++];
                next[c] = child_start[c];
                stack.push_back(c);
                continue;
            }
            stack.pop_back();
            post.push_back(k);
        }
    }
    for (Index k = 0; k < n; ++k)
        perm[k] = perm_[post[k]];
    perm_.swap(perm);
    for (Index k = 0; k < n; ++k)
        inv_perm_[perm_[k]] = k;
}
