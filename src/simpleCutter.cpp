#include "simpleCutter.h"

#include <memory>
#include <utility>

std::vector<VertexID> SimpleCutter::getGrowingCut(int cutoff) {
    // TODO handle disconnected case dist == -1.
    // Choose K root vertices and compute distances from them.
    if (!sampleRoots(K, 3, K)) // TODO handle mindist 2.
        return {};

    std::vector<VertexID> cut;
    // Reset inS, bounds, neigh, borders, subgraphs.
    for (VertexID v : graph.vertices())
        neigh[v] = 0;
    for (int i = 0; i < K; ++i) {
        for (VertexID v : graph.vertices())
            inS[i][v] = false;
        inS[i][roots[i]] = true;
        bounds[i].lower = bounds[i].upper = 1;
        boundedSize[i] = 1;
        borders[i].clear();
        assert(borders[i].capacity() == graph.nVertices() + 1);
        for (VertexID v : graph[roots[i]]) {
            assert(!neigh[v]);
            neigh[v] = i + 1;
            borders[i].push(v);
        }
        subgraphs[i].clear();
        subID[roots[i]] = subgraphs[i].addVertex();
    }

    while (true) {
        // Remove top border vertices that became cut vertices.
        for (int i = 0; i < K; ++i) {
            while (!borders[i].empty() && neigh[borders[i].front()] == -1)
                borders[i].pop();
            assert(borders[i].empty() || neigh[borders[i].front()] == i + 1);
        }
        // Pick a side, among sides with a non-empty border: minimizing bounds[side].upper.
        int side = -1;
        int upperBound = graph.nVertices() + 1;
        for (int i = 0; i < K; ++i) {
            if (borders[i].empty())
                continue;
            if (bounds[i].upper < upperBound) {
                upperBound = bounds[i].upper;
                side = i;
            }
        }
        // Return when there are no more sides with non-empty border.
        if (side == -1) {
            // log(-2, "Yield cut %zd td %d<=S<=%d \t %d<=T<=%d\n", cut.size(), bounds[0].lower, bounds[0].upper, bounds[1].lower, bounds[1].upper);
            return cut;
        }

        // Pick a pierce vertex p: first from the border queue (closest to root).
        VertexID p = borders[side].pop();
        assert(neigh[p] == side + 1);
        // Add p to S_side.
        inS[side][p] = true;
        ++bounds[side].upper;
        subID[p] = subgraphs[side].addVertex();
        for (VertexID v : graph[p])
            if (inS[side][v])
                subgraphs[side].addEdge(subID[p], subID[v]);
        // Update N(p) as border vertices, add to cut if already in another border.
        for (VertexID v : graph[p]) {
            if (!neigh[v]) {
                neigh[v] = side + 1;
                borders[side].push(v);
            } else if (neigh[v] != side + 1 && neigh[v] != -1) {
                cut.push_back(v);
                neigh[v] = -1;
                if (cut.size() >= cutoff)
                    return {};
            }
        }

        if (environment.signalled)
            return {};
    }
    __builtin_unreachable();
}


VertexID SimpleCutter::bfs(VertexID root, std::vector<int>& dist) {
    for (int& d : dist)
        d = -1;
    queue.clear();
    queue.push(root);
    dist[root] = 0;
    VertexID u = root;
    while (!queue.empty()) {
        u = queue.pop();
        int d = dist[u];
        for (const VertexID v : graph[u]) {
            if (dist[v] >= 0)
                continue;
            dist[v] = d + 1;
            queue.push(v);
        }
    }
    return u;
}


/** Find a cut by sampling two root vertices and returning the neighobrhood of vertices with dist[0][v] < dist[1][v] - 1. */
std::vector<VertexID> SimpleCutter::getDistanceCut(int cutoff) {
    std::vector<VertexID> cut;
    cut.reserve(cutoff);

    if (!sampleRoots(2, 2, 2))
        return {};

    for (VertexID u : graph.vertices()) {
        if (dists[0][u] == dists[1][u] || dists[0][u] == dists[1][u] - 1) {
            for (VertexID v: graph[u]) {
                if (dists[0][v] < dists[1][v] - 1) {
                    cut.push_back(u);
                    if (cut.size() >= cutoff)
                        return {};
                    break;
                }
            }
            if (environment.signalled)
                return {};
        }
    }

    return cut;
}


VertexID SimpleCutter::mostDistant(std::vector<VertexID> sources) {
    if (sources.empty())
        return getRandomVertex(graph, environment.randGen);
    queue.clear();
    std::vector<bool> visited(graph.nVertices() + 1, false);
    for (VertexID s : sources) {
        assert(s);
        queue.push(s);
        visited[s] = true;
    }
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


bool SimpleCutter::sampleRootsAttempt(int k, int minDistance, int pushes) {
    assert(k <= K);
    for (int i = 0; i < k; ++i) {
        if (i >= pushes)
            roots[i] = getRandomVertex(graph, environment.randGen);
        else {
            std::vector<VertexID> otherRoots;
            otherRoots.reserve(i);
            for (int j = 0; j < i; ++j)
                otherRoots.push_back(roots[j]);
            roots[i] = mostDistant(otherRoots);
        }
    }
    pushes -= k;
    for (int push = 0; push < pushes; ++push) {
        std::vector<VertexID> otherRoots;
        otherRoots.reserve(k - 1);
        int i = push % k;
        for (int j = 0; j < k; ++j)
            if (j != i)
                otherRoots.push_back(roots[j]);
        roots[i] = mostDistant(otherRoots);
    }

    for (int i = 0; i < k; ++i) {
        assert(roots[i]);
        for (int j = 0; j < i; ++j)
            if (dists[j][roots[i]] < minDistance)
                return false;
        bfs(roots[i], dists[i]);
    }
    return true;
}


bool SimpleCutter::sampleRoots(int k, int minDistance, int pushes, int maxRetries) {
    for (int r = 0; r < maxRetries && !environment.signalled; ++r)
        if (sampleRootsAttempt(k, minDistance, pushes))
            return true;
    return false;
}
