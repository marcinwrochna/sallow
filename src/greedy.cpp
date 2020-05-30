#include "greedy.h"

#include <algorithm>
#include <cassert>
#include <vector>
#include <valarray>
#include <gsl/gsl>

#include "minHeap.h"


int degeneracy(const BasicOneBasedGraph& graph) {
    std::vector<int> initVector;
    initVector.reserve(graph.nVertices());
    for (VertexID v : graph.vertices())
        initVector.push_back(graph[v].size());
    MinHeap heap(initVector);
    Ordering ordering;
    ordering.reserve(graph.nVertices());
    int maxDeg = 0;
    while (!heap.empty()) {
        int u, deg;
        std::tie(u, deg) = heap.pop();
        assert(!contains(ordering, u));
        ordering.push_back(u);
        if (deg > maxDeg)
            maxDeg = deg;
        for (const VertexID v : graph[u])
            heap.tryDecrementValue(v);
    }

    return maxDeg;
}


/** The "find" from the union-find data structure. */
static VertexID getRoot(VertexID v, std::vector<VertexID>& ancestor) {
    while (ancestor[v] != v) {
        assert(ancestor[v]);
        ancestor[v] = ancestor[ancestor[v]]; // Path halving, see <https://doi.org/10.1145%2F62.2160>.
        v = ancestor[v];
    }
    return v;
}


TreedepthResult superFastGreedy(const BasicOneBasedGraph& graph, const Ordering& ordering,
                           int lookahead, int cutoff,
                           Environment& environment) {
    int n = graph.nVertices();
    assert(ordering.size() == n);
    TreedepthResult result;
    result.ordering = ordering;
    result.parent = std::vector<VertexID>(n + 1, 0);
    result.heights = std::vector<VertexID>(n + 1, 1);
    // The algorithm is essentially the same as TreedepthResult::computeParents,
    // see there for details. Only here the new ordering is decided as we go.
    // We process vertices one by one and maintain a treedepth decomposition of
    // processed vertices in result.parent. To speed up mapping from a vertex to
    // the root in that tree, we maintain an ancestor, which will be updated as in
    // the union-find data structure.
    std::vector<VertexID> ancestor(n + 1, 0);
    // Unprocessed vertices have ancestor == 0.
    for (int i = n - 1; i >= 0; --i) {
        // Check the lookahead-many last non-processed vertices
        // and choose the one minimizing the max height of neighboring processed components.
        int bestJ = 0;
        int bestMaxHeight = n;
        for (int j = 0; j < std::min(lookahead, i + 1); ++j) {
            VertexID u = result.ordering[i - j];
            int maxHeight = 0;
            for (VertexID v : graph[u]) {
                if (ancestor[v]) {
                    int height = result.heights[getRoot(v, ancestor)];
                    assert(height);
                    if (height > maxHeight) {
                        maxHeight = height;
                        if (maxHeight + 1 >= cutoff || environment.signalled) {
                            result.clear();
                            return result;
                        }
                        if (maxHeight >= bestMaxHeight)
                            break;
                    }
                }
            }
            if (maxHeight < bestMaxHeight) {
                bestJ = j;
                bestMaxHeight = maxHeight;
            }
        }
        std::swap(result.ordering[i], result.ordering[i - bestJ]);
        VertexID u = result.ordering[i];
        ancestor[u] = u; // The newly processed vertex is becomes the new root.
        for (VertexID v : graph[u]) {
            if (ancestor[v]) {
                VertexID root = getRoot(v, ancestor);
                if (root != u) {
                    result.parent[root] = u;
                    ancestor[root] = u;
                }
            }
        }
        result.heights[u] = bestMaxHeight + 1;
    }
    assert(result.parent[result.ordering[0]] == 0);
    result.depth = result.heights[result.ordering[0]];
    return result;
}


/**
 * Merge-sort source into target, except for one element of source.
 * Assumes target and source are already sorted and source contains except.
 */
static void mergeUnique(std::vector<VertexID>& target, const std::vector<VertexID>& source, VertexID except) {
    std::vector<VertexID> newTarget;
    assert((long int)target.size() + (long int)source.size() - 1 >= 0);
    newTarget.reserve(target.size() + source.size() - 1);
    auto s = source.begin();
    auto t = target.begin();
    while (true) {
        if (s == source.end()) {
            for (; t != target.end(); ++t)
                newTarget.push_back(*t);
            break;
        } else if (t == target.end()) {
            for (; s != source.end(); ++s)
                if (*s != except)
                    newTarget.push_back(*s);
            break;
        } else if (*s < *t) {
            if (*s != except)
                newTarget.push_back(*s);
            ++s;
        } else if (*t < *s) {
            newTarget.push_back(*t);
            ++t;
        } else { // *s == *t
            newTarget.push_back(*s);
            ++t;
            ++s;
        }
    }
    target = std::move(newTarget);
}


/**
 * Merge-sort source into target, except for one element of source and one from target.
 * Assumes sExcept != tExcept, that they occur in source and target, respectively, and
 * that target and source are already sorted.
 * Returns the difference of resulting target.size() - original target.size().
 */
static int mergeUnique(std::vector<VertexID>& target, const std::vector<VertexID>& source,
                 VertexID sExcept, VertexID tExcept) {
    std::vector<VertexID> newTarget;
    assert((long int)target.size() + (long int)source.size() - 2 >= 0);
    newTarget.reserve(target.size() + source.size() - 2);
    auto s = source.begin();
    auto t = target.begin();
    while (true) {
       if (s == source.end()) {
            for (; t != target.end(); ++t)
                if (*t != tExcept)
                newTarget.push_back(*t);
            break;
        } else if (t == target.end()) {
            for (; s != source.end(); ++s)
                if (*s != sExcept)
                    newTarget.push_back(*s);
            break;
        } else if (*s < *t) {
            if (*s != sExcept)
                newTarget.push_back(*s);
            ++s;
        } else if (*t < *s) {
            if (*t != tExcept)
                newTarget.push_back(*t);
            ++t;
        } else { // *s == *t
            newTarget.push_back(*s);
            ++t;
            ++s;
        }
    }
    int result = (int)newTarget.size() - target.size();
    target = std::move(newTarget);
    return result;
}


/**
 * An implementation of greedy elimination which stores processed contracted components.
 *
 * This avoids storing dense graphs after elimination, but degrees of remaining vertices
 * can only be estimated. The lookup version turns out to work better.
 * We use this version when alpha=0, so degrees don't matter.
 */
static TreedepthResult fastGreedyElimByHeurVector(const BasicOneBasedGraph& graph,
                                 const std::vector<int>& centralities,
                                 int cutoff,
                                 int alpha, int beta, int gamma,
                                 Environment& environment) {
    int nVertices = graph.nVertices();
    TreedepthResult result;
    result.depth = 0;
    result.ordering.reserve(nVertices);
    result.parent = std::vector<VertexID>(nVertices + 1, 0);
    result.heights = std::vector<VertexID>(nVertices + 1, 1);
    // Make a working copy g of graph so that we can perform merges.
    std::vector<std::vector<int>> g(nVertices + 1);
    for (VertexID u : graph.vertices()) {
        g[u] = std::vector<int>(graph[u].begin(), graph[u].end());
        std::sort(g[u].begin(), g[u].end());
    }
    // Ancestors as in fastGreedy().
    std::vector<VertexID> ancestor(nVertices + 1, 0);

    std::vector<int> initVector;
    initVector.reserve(nVertices);
    for (VertexID u : graph.vertices())
        initVector.push_back(alpha * g[u].size() + beta * 1 + gamma * centralities[u]);
    MinHeap<int> heap(initVector, &environment.randGen);
    initVector.clear();

    while (!heap.empty()) {
        VertexID u = heap.pop().first;
        assert(!ancestor[u]);
        result.ordering.push_back(u);
        ancestor[u] = u; // The newly processed vertex is becomes the new root.
        for (VertexID v : graph[u]) {
            if (ancestor[v]) {
                VertexID root = getRoot(v, ancestor);
                if (root != u) {
                    result.parent[root] = u;
                    ancestor[root] = u;
                    assert(result.heights[u] >= result.heights[root] + 1);
                    // result.heights[u] = std::max(result.heights[u], result.heights[root] + 1);
                    mergeUnique(g[u], g[root], u, v);
                    g[root].clear(); // TODO check if other greedy should clear too?
                    g[root].shrink_to_fit(); // vector::clear always keeps its capacity.
                }
                else
                    mergeUnique(g[u], std::vector<int>({u}), u, v); // TODO think of something better.
            }
        }
        int lastV = 0;
        for (VertexID v : g[u]) {
            assert(!ancestor[v]);
            assert(v != lastV);
            lastV = v;
            result.heights[v] = std::max(result.heights[v], result.heights[u] + 1);
            int deg = std::max(g[v].size(), g[u].size()); // A heuristic, we can't afford merging all neighbors of v to g[v].
            heap.maximizeValue(v, alpha * deg + beta * result.heights[v] + gamma * centralities[v]);
        }
        if (result.heights[u] + g[u].size() >= cutoff || environment.signalled) {
            result.clear();
            return result;
        }
    }
    assert(result.ordering.size() == nVertices);
    std::reverse(result.ordering.begin(), result.ordering.end());
    assert(result.parent[result.ordering[0]] == 0);
    result.depth = result.heights[result.ordering[0]];
    return result;
}


/**
 * An implementation of greedy elimination which reevaluates degrees on the fly.
 *
 * Similarly as in fastGreedyElimByHeurVector, we store processed nodes as contracted components
 * instead of neighborhoods of the graph after elimination (which can be dense).
 * However, when before taking the next vertex from the minheap, we recompute
 * our evaluation of the top of the heap at most `lookahead` many times,
 * possibly pushing it down the heap.
 *
 * Lookahead 2 is a nice safely fast default.
 * Lookahead 1 is generally slower and worse, somehow.
 * Lookahead 1024 is slower but reasonable and almost as good as unbounded.
 */
static TreedepthResult lookGreedyElimByHeurVector(const BasicOneBasedGraph& graph,
                                 const std::vector<int>& centralities,
                                 int lookahead, int cutoff,
                                 int alpha, int beta, int gamma,
                                 Environment& environment) {
    int nVertices = graph.nVertices();
    TreedepthResult result;
    result.depth = 0;
    result.ordering.reserve(nVertices);
    result.parent = std::vector<VertexID>(nVertices + 1, 0);
    result.heights = std::vector<VertexID>(nVertices + 1, 1);
    // Make a working copy g of graph so that we can perform merges.
    std::vector<std::vector<int>> g(nVertices + 1);
    for (VertexID u : graph.vertices()) {
        g[u] = std::vector<int>(graph[u].begin(), graph[u].end());
        std::sort(g[u].begin(), g[u].end());
    }
    // Ancestors as in fastGreedy().
    std::vector<VertexID> ancestor(nVertices + 1, 0);

    std::vector<int> initVector;
    initVector.reserve(nVertices);
    for (VertexID u : graph.vertices())
        initVector.push_back(alpha * g[u].size() + beta * 1 + gamma * centralities[u]);
    MinHeap<int> heap(initVector, &environment.randGen);
    initVector.clear();

    while (!heap.empty()) {
        // Update heuristic evaluation of top of heap (at most lookahead-many times).
        VertexID u;
        std::vector<int> merged;
        for (int i = 0; true; ++i) {
            u = heap.top().first;
            assert(!ancestor[u]);
            merged = g[u]; // We do a copy that will often be discarded and recomputed; this is worth it when memory usage starts to slow us down.
            for (VertexID v : g[u]) {
                if (ancestor[v]) {
                    VertexID root = getRoot(v, ancestor); // MAYBE optimize if we have many neighbors in the same component (with the same root)? Take care of always removing v from merged, though.
                    // result.heights[u] = std::max(result.heights[u], result.heights[root] + 1);
                    assert(result.heights[u] >= result.heights[root] + 1);
                    mergeUnique(merged, g[root], u, v);
                }
            }
            if (i + 1 == lookahead)
                break;
            heap.maximizeValue(u, alpha * merged.size() + beta * result.heights[u] + gamma * centralities[u]);
            if (heap.top().first == u) // If top is still minimum after updating, it will always be.
                break;
        }

        heap.pop();
        result.ordering.push_back(u);
        ancestor[u] = u; // The newly processed vertex is becomes the new root.
        g[u] = std::move(merged);
        for (VertexID v : graph[u]) {
            if (ancestor[v]) {
                VertexID root = getRoot(v, ancestor); // MAYBE update getRoot path compression so that in this call we can always immediately return ancestor[v]?
                if (root != u) {
                    result.parent[root] = u;
                    ancestor[root] = u;
                    g[root].clear();
                    g[root].shrink_to_fit(); // vector::clear always keeps its capacity.
                }
            }
        }
        [[maybe_unused]] int lastV = 0;
        for (VertexID v : g[u]) {
            assert(!ancestor[v] && v != lastV);
            lastV = v;
            result.heights[v] = std::max(result.heights[v], result.heights[u] + 1);
            int deg = g[v].size() + g[u].size() - 1; // std::max(g[v].size(), g[u].size());
            heap.maximizeValue(v, alpha * deg + beta * result.heights[v] + gamma * centralities[v]);
        }
        if (result.heights[u] + g[u].size() >= cutoff || environment.signalled) {
            result.clear();
            return result;
        }
    }
    assert(result.ordering.size() == nVertices);
    std::reverse(result.ordering.begin(), result.ordering.end());
    assert(result.parent[result.ordering[0]] == 0);
    result.depth = result.heights[result.ordering[0]];
    return result;
}


/**
 * The standard implementation of greedy elimination.
 *
 * We compute the actual neighborhoods of vertices after each elimination.
 * This can make the graph pretty dense, greatly slowing down the algorithm
 * and possibly getting out of memory, if the expected treedepth is large.
 */
static TreedepthResult greedyElimByHeurVector(const BasicOneBasedGraph& graph,
                                 const std::vector<int>& centralities,
                                 int cutoff,
                                 int alpha, int beta, int gamma,
                                 Environment& environment) {
    int nVertices = graph.nVertices();
    TreedepthResult result;
    result.depth = 0;
    result.ordering.reserve(nVertices);
    // Make a working copy g of graph so that we can perform merges.
    std::vector<std::vector<int>> g(graph.nVertices() + 1);
    for (VertexID u : graph.vertices()) {
        g[u] = std::vector<int>(graph[u].begin(), graph[u].end());
        std::sort(g[u].begin(), g[u].end());
    }

    std::vector<int> initVector;
    initVector.reserve(nVertices);
    for (VertexID u : graph.vertices())
        initVector.push_back(alpha * g[u].size() + beta * 1 + gamma * centralities[u]);
    MinHeap<int> heap(initVector, &environment.randGen);
    std::vector<int> height(nVertices + 1, 1);
    int lastHeight = 0;
    while (!heap.empty()) {
        VertexID u = heap.pop().first;
        assert(!contains(result.ordering, u));
        result.ordering.push_back(u);
        lastHeight = height[u];
        if (lastHeight + g[u].size() >= cutoff || environment.signalled) {
            result.clear();
            return result;
        }
        for (VertexID v : g[u]) {
            assert(v != u);
            assert(!contains(result.ordering, v));
            assert(heap.contains(v));
            height[v] = std::max(height[v], lastHeight + 1);
            mergeUnique(g[v], g[u], v, u);
            heap.maximizeValue(v, alpha * g[v].size() + beta * height[v] + gamma * centralities[v]);
        }
        g[u].clear();
        g[u].shrink_to_fit();
    }
    assert(result.ordering.size() == nVertices);
    result.depth = lastHeight;
    std::reverse(result.ordering.begin(), result.ordering.end());
    result.computeParents(graph);
    return result;
}


TreedepthResult greedyElimByHeur(const BasicOneBasedGraph& graph,
                                 const std::vector<int>& centralities,
                                 int cutoff,
                                 int alpha, int beta, int gamma,
                                 Environment& environment) {
    assert(!gamma || centralities.size() == graph.nVertices() + 1);
    if (alpha == 0)
        return fastGreedyElimByHeurVector(graph, centralities, cutoff, alpha, beta, gamma, environment);
    bool useLookupVersion;
    if (graph.nVertices() >= 100'000)
        useLookupVersion = true;
    else if (cutoff >= 1000)
        useLookupVersion = true;
    else if (cutoff < 250)
        useLookupVersion = false;
    else
        useLookupVersion = environment.randGen() % 2;
    if (useLookupVersion)
        return lookGreedyElimByHeurVector(graph, centralities, 500, cutoff, alpha, beta, gamma, environment);
    else
        return greedyElimByHeurVector(graph, centralities, cutoff, alpha, beta, gamma, environment);
}


TreedepthResult locallyImprove(const BasicOneBasedGraph& graph, TreedepthResult& old, Environment& environment) {
    int nVertices = graph.nVertices();
    TreedepthResult result;
    result.depth = 0;
    result.ordering.reserve(nVertices);
    result.parent = std::vector<VertexID>(nVertices + 1, 0);
    result.heights = std::vector<VertexID>(nVertices + 1, 1);
    old.computeHeights();
    // A child v of a node p is critical if it is the only one with height[v] == height[p] - 1.
    std::vector<VertexID> critical(nVertices + 1, 0);
    for (VertexID u : graph.vertices())
        if (old.parent[u] && old.heights[u] == old.heights[old.parent[u]] - 1)
            critical[old.parent[u]] = u;
    for (VertexID u : graph.vertices())
        if (old.parent[u] && old.heights[u] == old.heights[old.parent[u]] - 1)
            if (critical[old.parent[u]] != u)
                critical[old.parent[u]] = 0;
    // Reorder so that parents come before children and the critical child of a node comes before its other children.
    Ordering ordering;
    ordering.reserve(graph.nVertices() + 1);
    for (VertexID u : old.ordering) {
        if (old.parent[u] && critical[old.parent[u]] == u) // If we're a critical child, we were already pushed.
            continue;
        for (VertexID v = u; v; v = critical[v]) // Push u, it's critical child, grandchild, and so on.
            ordering.push_back(v);
    }
    // Ancestors shortcutting result.parent as in superFastGreedy().
    std::vector<VertexID> ancestor(nVertices + 1, 0);

    for (int i = nVertices - 1; i >= 0; --i) {
        VertexID u = ordering[i];
        if (ancestor[u])
            continue; // u was already processed as the parent of its critical child.
        // If u is p's critical child, try swapping p with u.
        VertexID p = old.parent[u];
        if (p && critical[p] != u) // If we have a critical sibling, assert it was not processed yet.
            assert(!ancestor[critical[p]]);
        if (!p || critical[p] != u)
            goto noSwap;
        for (VertexID v : graph[p]) {
            if (!ancestor[v])
                continue;
            VertexID root = getRoot(v, ancestor);
            assert(root);
            // assert(old.parent[root] == p || old.parent[root] == u);
            assert(result.heights[root]);
            if (result.heights[root] + 2 >= old.heights[u] + 1)
                goto noSwap; // Don't swap if height with u and p would not improve.
        }
        // Process p.
        ancestor[p] = p;
        result.ordering.push_back(p);
        for (VertexID v : graph[p]) {
            if (!ancestor[v])
                continue; // Unprocessed vertices have ancestor == 0.
            VertexID root = getRoot(v, ancestor);
            if (root == p)
                continue;
            result.parent[root] = p;
            ancestor[root] = p;
            result.heights[p] = std::max(result.heights[p], result.heights[root] + 1);
        }
    noSwap:
        // Process u.
        ancestor[u] = u; // The newly processed vertex becomes the new root.
        result.ordering.push_back(u);
        for (VertexID v : graph[u]) {
            if (!ancestor[v])
                continue; // Unprocessed vertices have ancestor == 0.
            VertexID root = getRoot(v, ancestor);
            if (root == u)
                continue;
            result.parent[root] = u;
            ancestor[root] = u;
            result.heights[u] = std::max(result.heights[u], result.heights[root] + 1);
        }

        if (environment.signalled) {
            result.clear();
            return result;
        }
    }
    assert(result.ordering.size() == graph.nVertices());
    std::reverse(result.ordering.begin(), result.ordering.end());
    assert(result.parent[result.ordering[0]] == 0);
    result.depth = result.heights[result.ordering[0]];
    return result;
}


std::vector<VertexID> cutFromOrdering(const BasicOneBasedGraph& graph, const Ordering& ordering, int minSize) {
    int nVertices = graph.nVertices();
    assert(ordering.size() == nVertices);
    // The algorithm is essentially the same as TreedepthResult::computeParents,
    // we just keep track of sizes.
    std::vector<VertexID> ancestor(nVertices + 1, 0);
    std::vector<VertexID> size(nVertices + 1, 0);
    int nLargeComponents = 0;
    int bestIndex = 0;
    for (int i = nVertices - 1; i >= 0; --i) {
        VertexID u = ordering[i];
        ancestor[u] = u; // The newly processed vertex becomes the new root.
        size[u] = 1;
        if (size[u] >= minSize)
            ++nLargeComponents;
        for (VertexID v : graph[u]) {
            if (ancestor[v]) { // Unprocessed vertices have ancestor == 0.
                VertexID r = getRoot(v, ancestor);
                if (r != u) {
                    ancestor[r] = u;
                    if (size[r] >= minSize)
                        --nLargeComponents;
                    if (size[u] < minSize && size[u] + size[r] >= minSize)
                        ++nLargeComponents;
                    size[u] += size[r];
                }
            }
        }
        // After considering ordering[i,..,n-1], so like after removing ordering[0,...,i-1],
        // we have nLargeComponents components of size >= minSize.
        if (nLargeComponents >= 2)
            bestIndex = i;
    }
    assert(nLargeComponents == (nVertices >= minSize) ? 1 : 0);
    if (!bestIndex)
        return {};
    std::vector<VertexID> result;
    result.reserve(bestIndex);
    for (int i = 0; i < bestIndex; ++i)
        result.push_back(ordering[i]);
    return result;
}


TreedepthResult topDownGreedy(const BasicOneBasedGraph& graph, const Ordering& ordering) {
    int nVertices = graph.nVertices();
    std::vector<int> whichNode(nVertices + 1, 0);
    std::vector<std::vector<VertexID>> node(1); // node[0] is a dummy node at depth 0.
    node.reserve(nVertices / 2 + 1);
    std::vector<int> nodeDepth(1, 0);
    nodeDepth.reserve(nVertices / 2 + 1);
    // parent[node][0] is the parent, parent[node][logd] is the 2^logd-th ancestor.
    // So this is defined for logd < bitlength(nodeDepth).
    std::vector<std::vector<int>> parent(1);
    parent.reserve(nVertices / 2 + 1);
    for (int i = 0; i < nVertices; ++i) {
        VertexID u = ordering[i];
        // Order processed neighbors of u starting from the deepest.
        std::vector<std::pair<int, int>> neigh; // pairs <node,depth>.
        neigh.reserve(graph[u].size());
        for (VertexID v : graph[u])
            if (whichNode[v])
                neigh.emplace_back(whichNode[v], nodeDepth[whichNode[v]]);
        std::sort(neigh.begin(), neigh.end(), [](std::pair<int, int> a, std::pair<int, int> b){ return a.second > b.second; });
        int deepestNeigh = neigh.empty() ? 0 : neigh[0].second;
        // Check if all vertices in neigh are a chain of ancestors.
        bool isChain = true;
        for (int j = 0; j < (int)neigh.size() - 1; ++j) {
            assert(j + 1 < neigh.size());
            while (neigh[j].second > neigh[j + 1].second) {
                int d = neigh[j].second - neigh[j + 1].second;
                int logd = bitlength(d) - 1;
                d = 1 << logd;
                neigh[j].second -= d;
                assert(logd < parent[neigh[j].first].size());
                neigh[j].first = parent[neigh[j].first][logd];
            }
            if (neigh[j].first != neigh[j + 1].first) {
                isChain = false;
                break;
            }
        }
        if (isChain) {
            // Add a node j containing {u} below deepestNeigh.
            int j = node.size();
            whichNode[u] = j;
            node.emplace_back();
            node[j].push_back(u);
            assert(node[j].size() == 1);
            nodeDepth.push_back(nodeDepth[deepestNeigh] + 1);
            parent.emplace_back(bitlength(nodeDepth[deepestNeigh] + 1), 0);
            parent[j][0] = deepestNeigh;
            for (int logd = 0; (2 << logd) <= nodeDepth[deepestNeigh] + 1; ++logd) {
                assert(logd + 1 < parent[j].size());
                parent[j][logd + 1] = parent[parent[j][logd]][logd]; // TODO overflow here
            }
            continue;
        } else {
            // Compute LCA of all nodes.
            assert(!neigh.empty());
            int minDepth = neigh[neigh.size() - 1].second;
            for (int j = (int)neigh.size() - 2; j >= 0; --j) {
                // Make j the same depth as j+1 (which is minDepth).
                while (neigh[j].second > minDepth) {
                    int d = neigh[j].second - minDepth;
                    int logd = bitlength(d) - 1;
                    d = 1 << logd;
                    neigh[j].second -= d;
                    neigh[j].first = parent[neigh[j].first][logd];
                }
                // Move up until they're the same.
                while (neigh[j].first != neigh[j + 1].first) {
                    int logd = 0;
                    while ((2 << logd) < minDepth && parent[neigh[j].first][logd + 1] != parent[neigh[j + 1].first][logd + 1])
                        ++logd;
                    int d = 1 << logd;
                    assert(parent[neigh[j].first].size() == parent[neigh[j + 1].first].size());
                    assert(logd < parent[neigh[j].first].size());
                    neigh[j].first = parent[neigh[j].first][logd];
                    neigh[j + 1].first = parent[neigh[j + 1].first][logd];
                    minDepth -= d;
                }
                // Now neigh[j] is LCA of original neigh[j..], at depth minDepth.
            }
            whichNode[u] = neigh[0].first;
            node[neigh[0].first].push_back(u);
        }
    }

    TreedepthResult result;
    result.ordering.reserve(nVertices);
    for (int i = 0; i < node.size(); ++i)
        for (VertexID u : node[i])
            result.ordering.push_back(u);
    assert(result.ordering.size() == nVertices);
    result.computeParents(graph);
    result.computeHeights();
    result.depth = result.heights[result.ordering[0]];
    return result;
}


static int getSize(std::valarray<bool> a) {
    int result = 0;
    for (bool b : a)
        if (b)
            ++result;
    return result;
}


int greedyMinorBound(const BasicOneBasedGraph& graph, int prevLowerBound, Environment& environment) {
    int lowerBound = prevLowerBound;
    std::vector<std::valarray<bool>> g{graph.maxVertexID() + 1};
    for (VertexID u : graph.vertices())
        g[u] = std::valarray<bool>(false, graph.maxVertexID() + 1);
    int nEdges = 0;
    for (VertexID u : graph.vertices()) {
        for (VertexID v : graph[u]) {
            g[u][v] = true;
            ++nEdges;
        }
    }
    nEdges /= 2;
    std::vector<int> degree(graph.maxVertexID() + 1, 0);
    for (VertexID u : graph.vertices())
        degree[u] = graph[u].size();
    std::vector<int> initVector;
    initVector.reserve(graph.nVertices());
    for (VertexID u : graph.vertices())
        initVector.push_back(degree[u]);
    MinHeap<int> heap(initVector, &environment.randGen);
    std::valarray<bool> tmp(graph.maxVertexID() + 1);
    while (!heap.empty()) {
        auto top = heap.pop();
        VertexID u = top.first;
        int degU = top.second;
        lowerBound = std::max(lowerBound, degU + 1);
        while (2 * nEdges > static_cast<long long>(lowerBound - 1) * (2 * graph.nVertices() - lowerBound))
            ++lowerBound;
        assert(degU == degree[u]);
        // Let v be neighbour of u with smallest degree.
        VertexID v = 0;
        VertexID vDeg = graph.nVertices();
        int uDeg = 0;
        for (VertexID w : graph.vertices()) {
            if (g[u][w]) {
                if (degree[w] < vDeg) {
                    v = w;
                    vDeg = degree[v];
                }
                ++uDeg;
            }
        }
        assert(degree[u] == uDeg);
        if (!v)
            continue;
        // Contract u into v.
        tmp = g[u] && g[v];
        g[v] |= g[u];
        for (VertexID w : graph.vertices()) {
            if (g[u][w]) {
                g[w][u] = false;
                g[w][v] = true;
            }
        }
        g[v][v] = false;
        int intersection = 0;
        for (VertexID w : graph.vertices()) {
            if (tmp[w]) {
                heap.decreaseValue(w, --degree[w]);
                assert(degree[w] == getSize(g[w]));
                ++intersection;
            }
        }
        degree[v] += degree[u] - intersection - 2;
        heap.changeValueBy(v, degree[u] - intersection - 2);
        assert(degree[v] == getSize(g[v]));
        nEdges -= intersection + 1;
    }
    return lowerBound;
}
