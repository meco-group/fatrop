//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/elimination_order.hpp"

#include "fatrop/graph/block_sparsity.hpp"

#include <algorithm>
#include <limits>
#include <set>
#include <vector>

using namespace fatrop;

namespace
{
// Block-weighted minimum degree:
//   degree[v] = sum of sizes of v's current neighbours (excluding v itself).
// At each step we eliminate the vertex with smallest degree, merging all its
// neighbours into a clique. Ties are broken by lowest original index for
// determinism. Operates on a quotient elimination graph stored as std::set
// adjacency lists; this is the textbook implementation and is plenty fast
// for the modest block counts the solver typically sees. The cost is
// dominated by the clique creation step (O(deg^2) per elimination); when N
// is small this is negligible compared to the numerical factorisation.
struct EliminationGraph
{
    Index n;
    std::vector<std::set<Index>> adj;
    std::vector<Index> size;     // weight of each remaining super-vertex
    std::vector<Index> deg;      // cached weighted degree
    std::vector<bool> alive;

    explicit EliminationGraph(const BlockSparsityPattern &sp)
        : n(sp.num_blocks()), adj(sp.num_blocks()), size(sp.block_sizes()),
          deg(sp.num_blocks(), 0), alive(sp.num_blocks(), true)
    {
        for (Index v = 0; v < n; ++v)
        {
            const auto &nb = sp.neighbors(v);
            adj[v].insert(nb.begin(), nb.end());
        }
        for (Index v = 0; v < n; ++v)
            recompute_degree(v);
    }

    void recompute_degree(Index v)
    {
        Index d = 0;
        for (Index u : adj[v])
            d += size[u];
        deg[v] = d;
    }

    Index pick_min_degree() const
    {
        Index best = -1;
        Index best_deg = std::numeric_limits<Index>::max();
        for (Index v = 0; v < n; ++v)
        {
            if (!alive[v])
                continue;
            // Strict <: lowest original index wins ties (since we scan in order).
            if (deg[v] < best_deg)
            {
                best_deg = deg[v];
                best = v;
            }
        }
        return best;
    }

    void eliminate(Index v)
    {
        // Snapshot the neighbour set; we will modify adjacency.
        std::vector<Index> neigh(adj[v].begin(), adj[v].end());
        // Detach v from its neighbours.
        for (Index u : neigh)
            adj[u].erase(v);
        adj[v].clear();
        alive[v] = false;
        // Connect all pairs in neigh (clique fill-in).
        for (size_t a = 0; a < neigh.size(); ++a)
        {
            Index ua = neigh[a];
            for (size_t b = a + 1; b < neigh.size(); ++b)
            {
                Index ub = neigh[b];
                // std::set::insert is a no-op for existing edges.
                adj[ua].insert(ub);
                adj[ub].insert(ua);
            }
        }
        // Recompute degrees of touched vertices.
        for (Index u : neigh)
            recompute_degree(u);
    }
};
} // namespace

BlockEliminationOrder::BlockEliminationOrder(const BlockSparsityPattern &sp)
    : perm_(sp.num_blocks()), inv_perm_(sp.num_blocks())
{
    EliminationGraph g(sp);
    for (Index step = 0; step < sp.num_blocks(); ++step)
    {
        const Index v = g.pick_min_degree();
        perm_[step] = v;
        inv_perm_[v] = step;
        g.eliminate(v);
    }
}
