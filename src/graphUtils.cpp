#include "graphUtils.h"

#include <cassert>

#include "queue.h"


void TreedepthResult::clear() {
    depth = 0;
    ordering.clear();
    parent.clear();
    heights.clear();
    depths.clear();
}


void TreedepthResult::computeHeights() {
    if (!heights.empty())
        return;
    int n = ordering.size();
    heights = std::vector<int>(n + 1, 1);
    for (int i = n - 1; i > 0; --i) {
        VertexID u = ordering[i];
        heights[parent[u]] = std::max(heights[parent[u]], heights[u] + 1);
    }
}


void TreedepthResult::computeDepths() {
    if (!depths.empty())
        return;
    int n = ordering.size();
    depths = std::vector<int>(n + 1, 0);
    for (int i = 1; i < n; ++i) {
        VertexID u = ordering[i];
        depths[u] = depths[parent[u]] + 1;
    }
}


void TreedepthResult::computeParents(const BasicOneBasedGraph& graph) {
    if (!parent.empty())
        return;
    int n = graph.nVertices();
    assert(ordering.size() == n);
    parent = std::vector<VertexID>(n + 1, 0);
    // We process vertices from the end of `ordering`.
    // We maintain a treedepth decomposition of processed vertices in result.parent.
    // Thus each connected component of processed vertices has an associated tree.
    // To speed up mapping from a vertex to the root in that tree, we maintain an ancestor,
    // which will be updated as in the union-find data structure.
    // However, the root will be prescribed by us (not by size or rank), so the formal
    // guarantee is not Ackermann-style. Nevertheless, we expect those trees to have
    // small depth (at most `cutoff`): path compression in `ancestor` only improves that.
    std::vector<VertexID> ancestor(n + 1, 0);
    for (int i = n - 1; i >= 0; --i) {
        VertexID u = ordering[i];
        ancestor[u] = u; // The newly processed vertex becomes the new root.
        for (VertexID v : graph[u]) {
            if (ancestor[v]) { // Unprocessed vertices have ancestor == 0.
                VertexID r = v;
                while (ancestor[r] != r) {
                    assert(ancestor[r]);
                    ancestor[r] = ancestor[ancestor[r]]; // Path halving, see <https://doi.org/10.1145%2F62.2160>.
                    r = ancestor[r];
                }
                parent[r] = u;
                ancestor[r] = u;
            }
        }
    }
    parent[ordering[0]] = 0; // Could be == ordering[0] otherwise.
}


void TreedepthResult::assertCorrect(const BasicOneBasedGraph& graph) {
#ifndef NDEBUG
    if (!depth)
        return;
    assert(lowerBound <= depth);
    assert(parent.size() == graph.nVertices() + 1);
    assert(ordering.size() == graph.nVertices());
    // Check parents come before children in ordering.
    // (This implies `u → parent[u]` arrows define a forest oriented towards roots).
    std::vector<int> position(graph.nVertices() + 1);
    for (int i = 0; i < (int)ordering.size(); ++i)
        position[ordering[i]] = i;
    position[0] = -1;
    for ([[maybe_unused]] VertexID u : graph.vertices())
        assert(position[parent[u]] < position[u]);
    // Check there is exactly one root, namely ordering[0].
    for ([[maybe_unused]] VertexID u : graph.vertices())
        assert((!parent[u]) == (u == ordering[0]));

    // Check graph edges go between a vertex and its ancestor.
    for (VertexID u : graph.vertices()) {
        for (VertexID v : graph[u]) {
            if (position[u] > position[v])
                continue; // Only check u above v.
            bool ok = false;
            for (VertexID a = v; a != 0; a = parent[a]) {
                if (a == u) {
                    ok = true;
                    break;
                }
            }
            assert(ok);
        }
    }

    // Check depth of parent tree is exactly as declared.
    std::vector<int> tmpDepth(graph.nVertices() + 1, 99999999);
    tmpDepth[0] = 0; // Dummy has depth 0, root has depth 1.
    int maxDepth = 0;
    for (VertexID u : ordering) {
        tmpDepth[u] = tmpDepth[parent[u]] + 1;
        maxDepth = std::max(maxDepth, tmpDepth[u]);
    }
    assert(depth == maxDepth);
#endif
}


VertexID mostDistant(const BasicOneBasedGraph& graph, VertexID s) {
    Queue<VertexID> queue{graph.nVertices() + 1};
    std::vector<bool> visited(graph.maxVertexID() + 1, false);
    queue.push(s);
    visited[s] = true;
    VertexID u;
    while (!queue.empty()) {
        u = queue.pop();
        for (int v : graph[u]) {
            if (visited[v])
                continue;
            visited[v] = true;
            queue.push(v);
        }
    }
    return u;
}