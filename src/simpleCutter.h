#pragma once
#include <array>
#include <vector>

#include "environment.h"
#include "graph.h"
#include "graphUtils.h"
#include "queue.h"

class SimpleCutter
{
 public:
    constexpr static const int K = 2;

    SimpleCutter(const BasicOneBasedGraph& graph_, Environment& environment_)
        : graph(graph_),
          environment(environment_),
          neigh(std::vector<int>(graph.nVertices() + 1, 0)),
          subID(graph.nVertices() + 1, 0),
          queue(graph.nVertices() + 1)
    {
        for (int i = 0; i < K; ++i) {
            dists[i] = std::vector<int>(graph.nVertices() + 1, -1);
            inS[i] = std::vector<bool>(graph.nVertices() + 1, false);
            borders[i].reset(graph.nVertices() + 1);
        }
    }

    std::vector<VertexID> getGrowingCut(int cutoff);
    std::vector<VertexID> getDistanceCut(int cutoff);

 private:
    const BasicOneBasedGraph& graph;
    Environment& environment;

    std::array<VertexID, K> roots; // root[i] is the seed vertex from which S_i will grow.
    std::array<std::vector<int>, K> dists; // dist[i][v] is dist(root[i], v).
    std::array<std::vector<bool>, K> inS; // inS[i][v] is true iff v ∈ S_i.
    std::array<Bounds, K> bounds; // bounds on the treedepth of the subgraph induced by S_i.
    std::vector<int> neigh; // neigh[v] is -1 iff v ∈ cut, i+1 if v ∈ S_i ∪ N(S_i), 0 otherwise.
    std::array<Queue<VertexID>, K> borders; // border of each S_i, ordered by distance from root[i]
    std::array<BasicOneBasedGraph, K> subgraphs; // subgraphs idnuced by each S_i.
    std::vector<VertexID> subID; // subID[v] is the ID of v in the subgraph that contains it (or 0).
    std::array<int, K> boundedSize; // size of S_i when it was last bounded by bounds callback.

    Queue<int> queue; // Queue reused by all BFSs.

    VertexID bfs(VertexID root, std::vector<int>& dist);
    VertexID mostDistant(std::vector<VertexID> roots);
    bool sampleRoots(int k, int minDistance, int pushes, int maxRetries = 10);
    bool sampleRootsAttempt(int k, int minDistance, int pushes);
};
