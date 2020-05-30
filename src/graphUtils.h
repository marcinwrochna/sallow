#pragma once
#include <algorithm>
#include <vector>

#include "graph.h"

/**
 * A treedepth decomposition, with other useful stats and methods.
 * Depth == 0 means we have no meaningful result.
 * Otherwise at least parent and ordering are computed.
 * Heights and depths can be computed if needed, in almost linear time.
 * The root (that is, ordering[0]) has parent == 0.
 */
struct TreedepthResult
{
    int depth = 0;
    int lowerBound = 0;
    Ordering ordering; // Parents before children; root is ordering[0].
    std::vector<VertexID> parent;
    std::vector<int> heights; // Height of subtree rooted at given vertex; 1 for leaves.
    std::vector<int> depths;  // Length of path to root; 0 for root.

    /*
     * Set depth to 0 and clear all data.
     * Does not actually clear memory, vectors are still reserve to avoid reallocation.
     */
    void clear();

    /** Compute heights based on parents and ordering (if `heights` is empty). */
    void computeHeights();

    /** Compute depths based on parents and ordering (if `depths` is empty). */
    void computeDepths();

    /** Compute parents based on ordering (if `parent` is empty). */
    void computeParents(const BasicOneBasedGraph& graph);

    /**
     * Assert all of the decomposition is ok (depth, ordering, parent).
     *
     * This is a no-op if NDEBUG is defined.
     */
    void assertCorrect(const BasicOneBasedGraph& graph);
};



/** Return a vertex at maximum distance from s (via BFS).*/
VertexID mostDistant(const BasicOneBasedGraph& graph, VertexID s);

/** Return a vertex uniformly at random. */
template <typename RandGen>
VertexID getRandomVertex(const BasicOneBasedGraph& graph, RandGen& randGen) {
    std::uniform_int_distribution<VertexID> d{1, static_cast<VertexID>(graph.nVertices())};
    return d(randGen);
}

/**
 *  Return depth and vertices in dfs pre-order (parent before children).
 *  Only considers the connected component containing s.
 */
template <typename OneBasedGraph>
TreedepthResult dfsPreOrder(const OneBasedGraph& graph, int s) {
    TreedepthResult result;
    result.depth = 1;
    result.ordering.reserve(graph.nVertices());
    result.parent = std::vector<VertexID>(graph.maxVertexID() + 1, -1);
    int nVisited = 0;
    std::vector<int> depth(graph.maxVertexID() + 1, 0);
    typedef std::decay_t<decltype(graph[s].begin())> It;
    typedef std::decay_t<decltype(graph[s].end())> ItEnd;
    typedef std::tuple<VertexID, It, ItEnd> Node;
    std::vector<Node> stack(graph.nVertices());
    int se = 0;

    stack[se++] = std::make_tuple(s, graph[s].begin(), graph[s].end());
    depth[s] = 1;
    result.parent[s] = 0;
    result.ordering.push_back(s); // pre-order
    while (se) {
        VertexID u = std::get<0>(stack[se - 1]);
        It& b = std::get<1>(stack[se - 1]);
        ItEnd& e = std::get<2>(stack[se - 1]);
        if (b != e) {
            int v = *b;
            ++b;
            if (depth[v])
                continue;
            depth[v] = depth[u] + 1;
            result.parent[v] = u;
            result.depth = std::max(result.depth, depth[v]);
            result.ordering.push_back(v); // pre-order
            stack[se++] = std::make_tuple(v, graph[v].begin(), graph[v].end());
        } else {
            // Exiting node
            // result.ordering.push_back(u); // post-order
            --se;
            ++nVisited;
        }
    }
    // std::reverse(result.ordering.begin(), result.ordering.end()); // Post-order to elimination order.
    result.lowerBound = bitlength(result.depth); // The graph has a path of length depth.
    return result;
}
