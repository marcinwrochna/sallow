#pragma once
#include <cstdint>
#include <vector>

#include "graph.h"

struct TinyGraph
{
 public:
    int nVertices;
    std::uint64_t bitset;
    // Vertices are indexed 0..nVertices - 1.
    // Bitset records the presence of edges in order:
    // (0,1), (0,2), (1,2), (0,3), (1,3), (2,3), (0,4), (1,4), (2,4), (3,4), ...
    int getOffset(int u, int v) const noexcept {
        assert(0 <= u && u < nVertices);
        assert(0 <= v && v < nVertices);
        if (u > v)
            return u * (u - 1) / 2 + v;
        else
            return v * (v - 1) / 2 + u;
    }

    bool hasEdge(int u, int v) const noexcept {
        if (u == v)
            return false;
        return bitset & (1 << getOffset(u, v));
    };


    void addEdge(int u, int v) noexcept {
        assert(u != v);
        bitset = bitset | (1 << getOffset(u, v));
    };

    void removeIncidentEdges(int u) noexcept {
        for (int v = 0; v < nVertices; ++v)
            if (v != u)
                bitset = bitset & ~(1 << getOffset(u, v));
    }

    static TinyGraph getInducedSubgraph(const BasicOneBasedGraph& graph, const std::vector<int> vertices) noexcept {
        TinyGraph result{(int)vertices.size(), (std::uint64_t)0};
        std::vector<int8_t> newID(graph.nVertices() + 1, -1);
        for (int8_t u = 0; u < vertices.size(); ++u)
            newID[vertices[u]] = u;
        for (int8_t u = 0; u < vertices.size(); ++u)
            for (VertexID v : graph[vertices[u]])
                if (newID[v] != -1)
                    result.addEdge(u, newID[v]);
        return result;
    }

    static TinyGraph getInducedSubgraph(const TinyGraph& graph, const std::vector<int8_t> vertices) noexcept {
        TinyGraph result{(int)vertices.size(), (std::uint64_t)0};
        for (int8_t u = 0; u < vertices.size(); ++u)
            for (int8_t v = 0; v < u; ++v)
                if (graph.hasEdge(vertices[u], vertices[v]))
                    result.addEdge(u, v);
        return result;
    }
};


struct TinyTreedepthResult
{
    int depth;
    std::vector<int8_t> parent;
    std::vector<int8_t> ordering;
};


class TinyGraphDB
{
 public:
    #ifdef NDEBUG
        static const int maxSize = 7; // 2'097'152 graphs, about 150MiB total db size.
    #else
        static const int maxSize = 6;
    #endif

    const TinyTreedepthResult& getDecompositionConnected(const BasicOneBasedGraph& graph, std::vector<VertexID> vertices) const noexcept {
        assert(vertices.size() <= maxSize);
        TinyGraph tinyGraph = TinyGraph::getInducedSubgraph(graph, vertices);
        const TinyTreedepthResult& result = db[vertices.size()][tinyGraph.bitset];
        assert(result.depth >= 1 && result.depth <= vertices.size());
        return result;
    }

    void initialize();
 private:
    std::vector<TinyTreedepthResult> db[maxSize + 1];

    void initializeGraph(TinyGraph graph);
    void initializeConnectedGraph(TinyGraph graph);
};