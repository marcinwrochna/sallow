#include "tinyGraphDB.h"

#include <cstdio>

void TinyGraphDB::initialize() {
    fprintf(stderr, "Estimated DB size: %ld KiB\n", (sizeof(TinyTreedepthResult) + 3 * maxSize) * (1 << (maxSize * (maxSize - 1) / 2 - 10) ));
    db[1].resize(1);
    TinyTreedepthResult& t = db[1][0];
    t.depth = 1;
    t.parent.resize(1);
    t.parent[0] = -1;
    t.ordering.resize(1);
    t.ordering[0] = 0;
    for (int size = 2; size <= maxSize; ++size) {
        db[size].resize(1 << (size * (size - 1) / 2));
        for (std::uint64_t g = 0; g < 1 << (size * (size - 1) / 2); ++g)
            initializeGraph(TinyGraph{size, g});
    }
}


void TinyGraphDB::initializeGraph(TinyGraph graph) {
    TinyTreedepthResult& result = db[graph.nVertices][graph.bitset];
    result.parent.resize(graph.nVertices);
    result.ordering.reserve(graph.nVertices);
    result.depth = 0;
    std::vector<int8_t> stack;
    stack.reserve(graph.nVertices);
    std::vector<bool> visited(graph.nVertices);
    for (int8_t s = 0; s < graph.nVertices; ++s) {
        if (visited[s])
            continue;
        std::vector<int8_t> currentComponent;
        stack.push_back(s);
        visited[s] = true;
        while (!stack.empty()) {
            int8_t u = stack.back();
            stack.pop_back();
            currentComponent.push_back(u);
            for (int8_t v  = 0; v < graph.nVertices; ++v) {
                if (graph.hasEdge(u,v)) {
                    if (visited[v])
                    continue;
                    stack.push_back(v);
                    visited[v] = true;
                }
            }
        }
        if (currentComponent.size() == graph.nVertices)
            return initializeConnectedGraph(graph);
        TinyGraph subgraph = TinyGraph::getInducedSubgraph(graph, currentComponent);
        const TinyTreedepthResult& subResult = db[subgraph.nVertices][subgraph.bitset];
        assert(subResult.depth);
        result.depth = std::max(result.depth, subResult.depth);
        for (int v : subResult.ordering)
            result.ordering.push_back(currentComponent[v]);
        for (int v = 0; v < subgraph.nVertices; ++v)
            if (subResult.parent[v] != -1)
                result.parent[currentComponent[v]] = currentComponent[subResult.parent[v]];
            else
                result.parent[currentComponent[v]] = -1;
    }
}


void TinyGraphDB::initializeConnectedGraph(TinyGraph graph) {
    TinyTreedepthResult& result = db[graph.nVertices][graph.bitset];
    result.depth = graph.nVertices + 1;
    for (int8_t u = 0; u < graph.nVertices; ++u) {
        TinyGraph subgraph = graph;
        subgraph.removeIncidentEdges(u);
        const TinyTreedepthResult& r = db[subgraph.nVertices][subgraph.bitset];
        if (r.depth + 1 < result.depth) {
            result.depth = r.depth + 1;
            result.parent = r.parent;
            for (int8_t v = 0; v < graph.nVertices; ++v)
                if (result.parent[v] == -1)
                    result.parent[v] = u;
            result.parent[u] = -1;
            result.ordering.clear();
            result.ordering.push_back(u);
            for (int8_t v : r.ordering)
                if (v != u)
                    result.ordering.push_back(v);
        }
    }
}