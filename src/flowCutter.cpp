#include "flowCutter.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <limits>
#include <vector>

#include "graphUtils.h"
#include "utils.h"

// MAYBE replace vector with static containers (std::vector<int> is 6*sizeof(int), maybe gsl::dyn_array would work).

FlowCutter::FlowCutter(const BasicOneBasedGraph& graph_, Environment& environment_) :
    graph(graph_),
    environment(environment_),
    vertices(graph.nVertices() + 1),
    distS(graph.nVertices() + 1),
    distT(graph.nVertices() + 1),
    queue(graph.nVertices() + 1),
    cutoff(graph.nVertices()) {
    for (VertexID u : graph.vertices()) {
        vertices[u].neigh.reserve(graph[u].size());
        for (VertexID v : graph[u]) {
            assert(v >= 1 && v <= graph.nVertices());
            vertices[u].neigh.push_back(v);
        }
    }
}


FlowCutter::EndReason FlowCutter::findCuts(int _s, int _t, int _cutoff, YieldCutCallback yc) {
    s = _s;
    t = _t;
    cutoff = _cutoff;
    yieldCutCallback = yc;
    if (s == t)
        return EndReason::badST;
    for (const int v : vertices[s].neigh)
        if (v == t)
            return EndReason::badST;
    for (FlowCutter::Vertex& vertex : vertices) {
        vertex.inS = vertex.inT = false;
        vertex.boundaryS = vertex.boundaryT = false;
        vertex.inflow = vertex.outflow = 0;
        vertex.parentS = vertex.parentT = 0;
    }
    initDistS();
    initDistT();
    sizeS = sizeT = 0;
    clearedParentS = clearedParentT = true;
    vertices[s].parentS = s;
    addToS(s);
    vertices[t].parentT = t;
    addToT(t);
    flowEndpointsS.clear();
    flowEndpointsT.clear();
    boundsS.lower = boundsS.upper = boundsT.lower = boundsT.upper = 1;
    sizeBoundedS = sizeBoundedT = 1;
    boundsReachS.lower = boundsReachS.upper = boundsReachT.lower = boundsReachT.upper = 1;
    sizeBoundedReachS = sizeBoundedReachT = 1;

    // MAYBE try starting with a larger greedy flow.

    bool pierceSideIsS = true;
    pierceS = s;
    pierceT = t;
    // Invariant: everything reachable from S (in residual digraph) is already in S or reachable from pierceS.
    // In other words, neighbors of S are either in S, on flow-paths, or neighbors of pierceS. Same for T.
    // Neighbors of pierceS/T may postpone being be marked as boundaryS/T.
    while (true) {
        // Augment the flow until we have a maximum flow.
        if (pierceSideIsS) {
            int v;
            while ((v = bfsGrowParentS())) {
                if (flowEndpointsS.size() + 1 + boundsS.lower >= cutoff)
                    return EndReason::cutoff;
                if (environment.signalled)
                    return EndReason::signal;
                augmentFlowViaS(v);
            }
            v = bfsGrowParentT();
            assert(!v);
        } else { // Symmetric case.
            int v;
            while ((v = bfsGrowParentT())) {
                if (flowEndpointsT.size() + 1 + boundsT.lower >= cutoff)
                    return EndReason::cutoff;
                if (environment.signalled)
                    return EndReason::signal;
                augmentFlowViaT(v);
            }
            v = bfsGrowParentS();
            assert(!v);
        }
        // Push the cuts towards each other until we can't avoid increasing the flow.
        while (true) {
            // Now the flow is maximum; parentS, parentT are valid BFS structures from S,T.
            if (environment.signalled)
                return EndReason::signal;
            #ifndef NDEBUG
                checkInvariants();
            #endif
            pierceSideIsS = choosePierceSide();
            if (pierceSideIsS) {
                bfsGrowInS();
                // Now S == reachable from S; since the flow is maximum, flowEndpointsS is a cut.
                yieldSoftCut(flowEndpointsS);
                EndReason r = selectNewPierceS();
                if (r != EndReason::none) {
                    return r;
                } else if (vertices[pierceS].parentT) {
                    augmentFlowViaT(-vertices[pierceS].parentT);
                    break;
                } else { // Flow is still maximum. This is the most frequent branch.
                    [[maybe_unused]] int v = bfsGrowParentS();
                    assert(!v);
                }
            } else { // Symmetric case.
                bfsGrowInT();
                yieldSoftCut(flowEndpointsT);
                EndReason r = selectNewPierceT();
                if (r != EndReason::none) {
                    return r;
                } else if (vertices[pierceT].parentS) {
                    augmentFlowViaS(-vertices[pierceT].parentS);
                    break;
                } else {
                    [[maybe_unused]] int v = bfsGrowParentT();
                    assert(!v);
                }
            }
        }
    };
    assert(false);
}


bool FlowCutter::choosePierceSide() {
    int depthReachS = sizeReachS - sizeBoundedReachS + boundsReachS.upper;
    int depthReachT = sizeReachT - sizeBoundedReachT +  boundsReachT.upper;
    return depthReachS <= depthReachT;
}


void FlowCutter::yieldSoftCut([[maybe_unused]] std::vector<int> cut) {
    //printf("cut = "); for (int x : cut) printf(" %d", x); printf("\n");
    assert(checkIsCut(cut, s, t));
}


FlowCutter::EndReason FlowCutter::yieldCutS(EndReason r) {
    assert(checkIsCut(flowEndpointsS, s, t));
    if (r == EndReason::none && flowEndpointsS.size() + 1 + boundsS.lower >= cutoff)
        r = EndReason::cutoff;
    boundsS.upper += sizeS - sizeBoundedS;
    sizeBoundedS = sizeS;
    return yieldCutCallback(flowEndpointsS, r, sizeS);
}


FlowCutter::EndReason FlowCutter::yieldCutT(EndReason r) {
    assert(checkIsCut(flowEndpointsT, s, t));
    if (r == EndReason::none && flowEndpointsT.size() + 1 + boundsT.lower >= cutoff)
        r = EndReason::cutoff;
    boundsT.upper += sizeT - sizeBoundedT;
    sizeBoundedT = sizeT;
    return yieldCutCallback(flowEndpointsT, r, sizeT);
}



FlowCutter::EndReason FlowCutter::selectNewPierceS() {
    int result = -1;
    int bestDist = std::numeric_limits<int>::min();
    for (int i = 0; i < flowEndpointsS.size(); ++i) {
        int u = flowEndpointsS[i];
        if (!vertices[u].parentT && !vertices[u].boundaryT) {
            result = i;
            break;
        }
        // Neighbors of pierceT in boundaryS have parentT == -pierceT.
        if (vertices[u].boundaryT || vertices[u].parentT == -pierceT)
            continue;
        assert(vertices[u].parentT < 0);
        assert(!contains(vertices[u].neigh, pierceT));
        int dist = distT[u] - distS[u]; // The piercing heuristic.
        if (dist > bestDist) {
            bestDist = dist;
            result = i;
        }
    }
    if (result == -1) { // All flowEndpointsS are neighbors of T, no more cuts.
        pierceS = 0;
        return yieldCutS(EndReason::complete);
    }
    pierceS = flowEndpointsS[result];
    EndReason r = EndReason::none;
    if (vertices[pierceS].parentT) // The flow will be augmented, yield current cut before augmenting.
        r = yieldCutS(r);
    addToS(pierceS);
    pushFlowEndpointS(result);
    return r;
}


FlowCutter::EndReason FlowCutter::selectNewPierceT() {
    int result = -1;
    int bestDist = std::numeric_limits<int>::min();
    for (int i = 0; i < flowEndpointsT.size(); ++i) {
        int u = flowEndpointsT[i];
        if (!vertices[u].parentS && !vertices[u].boundaryS) {
            result = i;
            break;
        }
        // Neighbors of pierceS in boundaryT have parentS == -pierceS.
        if (vertices[u].boundaryS || vertices[u].parentS == -pierceS)
            continue;
        assert(vertices[u].parentS < 0);
        assert(!contains(vertices[u].neigh, pierceS));
        int dist = distS[u] - distT[u]; // The piercing heuristic.
        if (dist > bestDist) {
            bestDist = dist;
            result = i;
        }
    }
    if (result == -1) {
        pierceT = 0;
        return yieldCutT(EndReason::complete);
    }
    pierceT = flowEndpointsT[result];
    EndReason r = EndReason::none;
    if (vertices[pierceT].parentS)
        r = yieldCutT(r);
    addToT(pierceT);
    pushFlowEndpointT(result);
    return r;
}


/** Add vertex to S and update everything. */
void FlowCutter::addToS(int v) {
    assert(!vertices[v].inS);
    assert(vertices[v].boundaryS || vertices[v].parentS);
    sizeS++;
    vertices[v].inS = true;
    vertices[v].boundaryS = true;
    vertices[v].parentS = v;
}


/** Add vertex to T and update everything. */
void FlowCutter::addToT(int v) {
    assert(!vertices[v].inT);
    assert(vertices[v].boundaryT || vertices[v].parentT);
    sizeT++;
    vertices[v].inT = true;
    vertices[v].boundaryT = true;
    vertices[v].parentT = v;
}


void FlowCutter::pushFlowEndpointS(int index) {
    int oldEnd = flowEndpointsS[index];
    int newEnd = vertices[oldEnd].outflow;
    assert(newEnd && newEnd != oldEnd);
    vertices[oldEnd].inflow = vertices[oldEnd].outflow = 0;
    flowEndpointsS[index] = newEnd;
    vertices[newEnd].inflow = newEnd;
    vertices[newEnd].boundaryS = true;
}


void FlowCutter::pushFlowEndpointT(int index) {
    int oldEnd = flowEndpointsT[index];
    int newEnd = vertices[oldEnd].inflow;
    assert(newEnd && newEnd != oldEnd);
    vertices[oldEnd].outflow = vertices[oldEnd].inflow = 0;
    flowEndpointsT[index] = newEnd;
    vertices[newEnd].outflow = newEnd;
    vertices[newEnd].boundaryT = true;
}


/** Grow BFS from pierceS. Return boundary of T if found (giving an ST-path), zero otherwise. */
int FlowCutter::bfsGrowParentS() {
    clearedParentS = false;
    assert(!vertices[pierceS].parentT);
    assert(vertices[pierceS].inS);

    sizeReachS = sizeS - 1;  // Count pierceS only once, when popping from queue.
    boundsReachS = boundsS; // Reset bounds on reachS.
    sizeBoundedReachS = sizeBoundedS;
    queue.clear();
    queue.push(pierceS);
    while (!queue.empty()) {
        int u = queue.pop();
        ++sizeReachS;
        assert(vertices[u].parentS);
        assert(!vertices[u].parentT);
        for (const int v : vertices[u].neigh) {
            if (vertices[v].parentS)
                continue;
            assert(v == pierceT || !vertices[v].inT);

            // If u->v enters a flow-path, the residual digraph has an edge from u to v's inflow w instead.
            int w = vertices[v].inflow;
            if (!w) {
                if (v == pierceT)
                    return u;
                vertices[v].parentS = u;
                queue.push(v);
            } else {
                assert(v != pierceT);
                // When we enter a flow-path, we visit it all.
                assert(vertices[u].inflow != v); // If v->u is in flow, we already visited v.
                vertices[v].parentS = -u; // The intermediate vertex v is pre-reachable.
                if (w == u || w == v) // If u->v is in flow or v starts a flow-path, we skip it.
                    continue;
                assert(!vertices[w].parentT); // Otherwise already vertices[v].parentT > 0.
                int p = u;
                while (!vertices[w].parentS && vertices[w].inflow != w) {
                    vertices[w].parentS = p;
                    queue.push(w);
                    p = w;
                    w = vertices[w].inflow;
                }
                assert(!(vertices[w].parentS > 0));
                vertices[w].parentS = p;
                queue.push(w);
            }
        }
    }
    return 0;
}


void FlowCutter::bfsGrowInS() {
    queue.clear();
    queue.push(pierceS);
    while (!queue.empty()) {
        int u = queue.pop();
        for (const int v : vertices[u].neigh) {
            if (vertices[v].boundaryS)
                continue;
            // If u->v enters a flow-path, the residual digraph has an edge from u to v's inflow w instead.
            int w = vertices[v].inflow;
            if (!w) {
                addToS(v);
                queue.push(v);
            } else {
                assert(vertices[u].inflow != v);
                vertices[v].boundaryS = true;
                if (w == u || w == v)
                    continue;
                while (!vertices[w].boundaryS && vertices[w].inflow != w) {
                    addToS(w);
                    queue.push(w);
                    w = vertices[w].inflow;
                }
                assert(!vertices[w].inS);
                addToS(w);
                queue.push(w);
            }
        }
    }

    // Push flowEndpointsS.
    for (int& u : flowEndpointsS) {
        while (vertices[u].inS) {
            int old = u;
            u = vertices[u].outflow;
            vertices[old].inflow = vertices[old].outflow = 0;
        }
        assert(vertices[u].boundaryS);
        vertices[u].inflow = u;
    }

    // Now S == reachS, so we can update bounds.
    if (sizeBoundedReachS >= sizeBoundedS) {
        boundsS.upper += sizeBoundedReachS - sizeBoundedS;
        sizeBoundedS = sizeBoundedReachS;
        boundsS.updateUpper(boundsReachS.upper);
        boundsS.updateLower(boundsReachS.lower);
        boundsReachS = boundsS;
    }
}


/** Symmetric version of `bfsGrowParentS`, swapping S with T, inflow with outflow. */
int FlowCutter::bfsGrowParentT() {
    clearedParentT = false;
    assert(!vertices[pierceT].parentS);
    assert(vertices[pierceT].inT);

    sizeReachT = sizeT - 1;  // Count pierceT only when popping from queue.
    boundsReachT = boundsT; // Reset bounds on reachT.
    sizeBoundedReachT = sizeBoundedT;
    queue.clear();
    queue.push(pierceT);
    while (!queue.empty()) {
        int u = queue.pop();
        ++sizeReachT;
        assert(vertices[u].parentT);
        assert(!vertices[u].parentS);
        for (const int v : vertices[u].neigh) {
            if (vertices[v].parentT)
                continue;
            assert(v == pierceS || !vertices[v].inS);
            // If u->v enters a flow-path, the residual digraph has an edge from u to v's outflow w instead.
            int w = vertices[v].outflow;
            if (!w) {
                if (v == pierceS)
                    return u;
                vertices[v].parentT = u;
                queue.push(v);
            } else {
                assert(v != pierceS);
                // When we enter a flow-path, we visit it all.
                assert(vertices[u].outflow != v); // If u->v is in flow, we already visited v.
                vertices[v].parentT = -u;
                if (w == u || w == v) // If v->u is in flow or v ends a flow-path, we skip it.
                    continue;
                assert(!vertices[w].parentS); // Otherwise already vertices[v].parentS > 0.
                int p = u;
                while (!vertices[w].parentT && vertices[w].outflow != w) {
                    vertices[w].parentT = p;
                    queue.push(w);
                    p = w;
                    w = vertices[w].outflow;
                }
                assert(!(vertices[w].parentT > 0));
                vertices[w].parentT = p;
                queue.push(w);
            }
        }
    }
    return 0;
}


/** Symmetric version of `bfsGrowInS`, swapping S with T, inflow with outflow. */
void FlowCutter::bfsGrowInT() {
    queue.clear();
    queue.push(pierceT);
    while (!queue.empty()) {
        int u = queue.pop();
        for (const int v : vertices[u].neigh) {
            if (vertices[v].boundaryT)
                continue;
            // If u->v enters a flow-path, the residual digraph has an edge from u to v's outflow w instead.
            int w = vertices[v].outflow;
            if (!w) {
                addToT(v);
                queue.push(v);
            } else {
                assert(vertices[u].outflow != v);
                vertices[v].boundaryT = true;
                if (w == u || w == v)
                    continue;
                while (!vertices[w].boundaryT && vertices[w].outflow != w) {
                    addToT(w);
                    queue.push(w);
                    w = vertices[w].outflow;
                }
                assert(!vertices[w].inT);
                addToT(w);
                queue.push(w);
            }
        }
    }

    // Push flowEndpointsT.
    for (int& u : flowEndpointsT) {
        while (vertices[u].inT) {
            int old = u;
            u = vertices[u].inflow;
            vertices[old].inflow = vertices[old].outflow = 0;
        }
        assert(vertices[u].boundaryT);
        vertices[u].outflow = u;
    }

    // Now T == reachT, so we can update bounds.
    if (sizeBoundedReachT >= sizeBoundedT) {
        boundsT.upper += sizeBoundedReachT - sizeBoundedT;
        sizeBoundedT = sizeBoundedReachT;
        boundsT.updateUpper(boundsReachT.upper);
        boundsT.updateLower(boundsReachT.lower);
        boundsReachT = boundsT;
    }
}


/** Augment the flow given a neighbor of pierceT returned by `bfsGrowParentS`. Clears parentS/T. */
void FlowCutter::augmentFlowViaS(int last) {
    assert(last > 0);
    assert(!vertices[last].inS && !vertices[last].inT);
    // Follow parentS of last until we're in S.
    int a = last;
    int outflow = vertices[a].outflow; // The original outflow of a (before augmenting).
    vertices[a].outflow = a;
    for (int b = vertices[a].parentS; b != pierceS; a = b, b = vertices[b].parentS) {
        assert(!vertices[b].inS);
        assert(b > 0);
        if (!outflow) {
            outflow = vertices[b].outflow;
            vertices[a].inflow = b;
            vertices[b].outflow = a;
        } else {
            assert((vertices[b].inflow == a) == (outflow == b));
            while (outflow == b) {
                outflow = vertices[b].outflow;
                vertices[b].inflow = 0;
                vertices[b].outflow = 0;
                a = b;
                b = vertices[a].parentS;
                assert(b > 0);
            }
            // We're exiting an original flow-path.
            // Since we're going backwards, the BFS entered the path from b at tmp
            // and the residual digraph had an arc directly to its inflow a.
            if (b == pierceS)
                break;
            int tmp = outflow;
            assert(vertices[tmp].inflow == a);
            outflow = vertices[b].outflow;
            vertices[tmp].inflow = b;
            vertices[b].outflow = tmp;
        }
    }
    int first = a;
    if (outflow)
        first = outflow;
    vertices[first].inflow = first;
    flowEndpointsS.push_back(first);
    flowEndpointsT.push_back(last);
    assert(!vertices[first].inS);
    assert(!vertices[last].inT);
    assert(vertices[first].boundaryS == false);
    assert(vertices[last].boundaryT == false);
    vertices[first].boundaryS = true;
    vertices[last].boundaryT = true;
    clearParentS();
    clearParentT();
}


/** Augment the flow given a neighbor of pierceS returned by `bfsGrowParentT`. Clears parentS/T. */
void FlowCutter::augmentFlowViaT(int first) {
    assert(first > 0);
    assert(!vertices[first].inS && !vertices[first].inT);
    // Follow parentT of first until we're in T.
    int a = first;
    int inflow = vertices[a].inflow; // The original inflow of a (before augmenting).
    vertices[a].inflow = a;
    for (int b = vertices[a].parentT; b != pierceT; a = b, b = vertices[b].parentT) {
        assert(!vertices[b].inT);
        assert(b > 0);
        if (!inflow) {
            inflow = vertices[b].inflow;
            vertices[a].outflow = b;
            vertices[b].inflow = a;
        } else {
            assert((vertices[b].outflow == a) == (inflow == b));
            while (inflow == b) {
                inflow = vertices[b].inflow;
                vertices[b].outflow = 0;
                vertices[b].inflow = 0;
                a = b;
                b = vertices[a].parentT;
                assert(b > 0);
            }
            if (b == pierceT)
                break;
            int tmp = inflow;
            assert(vertices[tmp].outflow == a);
            inflow = vertices[b].inflow;
            vertices[tmp].outflow = b;
            vertices[b].inflow = tmp;
        }
    }
    int last = a;
    if (inflow)
        last = inflow;
    vertices[last].outflow = last;
    flowEndpointsS.push_back(first);
    flowEndpointsT.push_back(last);
    assert(!vertices[first].inS);
    assert(!vertices[last].inT);
    assert(vertices[first].boundaryS == false);
    assert(vertices[last].boundaryT == false);
    vertices[first].boundaryS = true;
    vertices[last].boundaryT = true;
    clearParentS();
    clearParentT();
}


void FlowCutter::clearParentS() {
    if (clearedParentS)
        return;
    clearedParentS = true;
    for (int u = 1; u < vertices.size(); ++u)
        if (!vertices[u].inS)
            vertices[u].parentS = 0;
}


void FlowCutter::clearParentT() {
    if (clearedParentT)
        return;
    clearedParentT = true;
    for (int u = 1; u < vertices.size(); ++u)
        if (!vertices[u].inT)
            vertices[u].parentT = 0;
}


void FlowCutter::initDistS() {
    for (int& d : distS)
        d = -1;
    queue.clear();
    queue.push(s);
    distS[s] = 0;
    while (!queue.empty()) {
        int u = queue.pop();
        int dist = distS[u];
        while (dist >= layerSizeS.size())
            layerSizeS.push_back(0);
        ++layerSizeS[dist];
        for (const int v : vertices[u].neigh) {
            if (distS[v] >= 0)
                continue;
            distS[v] = dist + 1;
            queue.push(v);
        }
    }
    // fprintf(stderr, "layerSizeS: %s\n", toString(layerSizeS).c_str());
}


void FlowCutter::initDistT() {
    for (int& d : distT)
        d = -1;
    queue.clear();
    queue.push(t);
    distT[t] = 0;
    while (!queue.empty()) {
        int u = queue.pop();
        int dist = distT[u];
        while (dist >= layerSizeT.size())
            layerSizeT.push_back(0);
        ++layerSizeT[dist];
        for (const int v : vertices[u].neigh) {
            if (distT[v] >= 0)
                continue;
            distT[v] = dist + 1;
            queue.push(v);
        }
    }
    // fprintf(stderr, "layerSizeT: %s\n", toString(layerSizeT).c_str());
}


void FlowCutter::checkInvariants() {
    for (int u = 1; u < vertices.size(); ++u) {
        const Vertex& U = vertices[u];
        if (u == pierceS)
            assert(U.inS);
        if (u == pierceT)
            assert(U.inT);
        // S,T are disjoint, boundaryS/T is exactly S/T or neighbors of S-pierceS/T-pierceT.
        if (U.inS) {
            assert(!U.inT);
            assert(!U.boundaryT);
            assert(!U.parentT);
            assert(U.parentS == u);
            if (u != pierceS)
                for ([[maybe_unused]] const int v : U.neigh)
                    assert(vertices[v].boundaryS);
            assert(!U.inflow);
            assert(!U.outflow);
        } else if (U.boundaryS) {
            assert(!U.inT);
            assert(!(U.parentT > 0));
            [[maybe_unused]] bool ok = false;
            for (const int v : U.neigh)
                if (vertices[v].inS)
                    ok = true;
            assert(ok);
            assert(U.inflow == u);
            assert(U.outflow);
        }
        if (U.inT) {
            assert(!U.inS);
            assert(!U.boundaryS);
            assert(!U.parentS);
            assert(U.parentT == u);
            if (u != pierceT)
                for ([[maybe_unused]] const int v : U.neigh)
                    assert(vertices[v].boundaryT);
            assert(!U.inflow);
            assert(!U.outflow);
        } else if (U.boundaryT) {
            assert(!U.inS);
            assert(!(U.parentS > 0));
            [[maybe_unused]] bool ok = false;
            for (const int v : U.neigh)
                if (vertices[v].inT)
                    ok = true;
            assert(ok);
            assert(U.outflow == u);
            assert(U.inflow);
        }

        // parentS, parentT are valid BFS trees from pierceS/pierceT
        if (!U.boundaryS) {
            int v = u;
            int w = vertices[v].parentS;
            if (w) {
                if (w < 0) {
                    w = -w;
                    assert(vertices[v].inflow && vertices[v].outflow);
                    assert(vertices[w].inflow != v && vertices[v].outflow != w);
                    assert(vertices[v].inflow != w && vertices[w].outflow != v);
                }
                while (!vertices[w].inS) {
                    v = w;
                    w = vertices[w].parentS;
                    assert(w > 0);
                }
                assert(w == pierceS);
            }
        }
        if (!U.boundaryT) {
            int v = u;
            int w = vertices[v].parentT;
            if (w) {
                if (w < 0) {
                    w = -w;
                    assert(vertices[v].outflow && vertices[v].inflow);
                    assert(vertices[w].outflow != v && vertices[v].inflow != w);
                    assert(vertices[v].outflow != w && vertices[w].inflow != v);
                }
                while (!vertices[w].inT) {
                    v = w;
                    w = vertices[w].parentT;
                    assert(w > 0);
                }
                assert(w == pierceT);
            }
        }

        // Flow-paths inflow and outflow describe consistent paths.
        if (!U.inS && !U.inT) {
            if (!U.inflow) {
                assert(!U.outflow);
            } else {
                if (U.inflow == u) {
                    [[maybe_unused]] bool ok = false;
                    for (const int v: flowEndpointsS)
                        if (u == v)
                            ok = true;
                    assert(ok);
                } else
                    assert(vertices[U.inflow].outflow == u);
                if (U.outflow == u) {
                    [[maybe_unused]] bool ok = false;
                    for (const int v: flowEndpointsT)
                        if (u == v)
                            ok = true;
                    assert(ok);
                } else
                    assert(vertices[U.outflow].inflow == u);
                // Check for cycles.
                int v = u;
                while (vertices[v].inflow != v)
                    v = vertices[v].inflow;
                v = u;
                while (vertices[v].outflow != v)
                    v = vertices[v].outflow;
            }
        }
    }
    assert(flowEndpointsS.size() == flowEndpointsT.size());
    for ([[maybe_unused]] const int u: flowEndpointsS) {
        assert(vertices[u].inflow == u);
        assert(!vertices[u].inS);
        assert(vertices[u].boundaryS);
    }
    for ([[maybe_unused]] const int u: flowEndpointsT) {
        assert(vertices[u].outflow == u);
        assert(!vertices[u].inT);
        assert(vertices[u].boundaryT);
    }
    // Sizes are maintained correctly.
    int realSizeS = 0;
    int realSizeT = 0;
    int realSizeReachS = 0;
    int realSizeReachT = 0;
    for (int u = 1; u < vertices.size(); ++u) {
        const Vertex& U = vertices[u];
        if (U.inS)
            realSizeS++;
        if (U.inT)
            realSizeT++;
        if (U.parentS > 0)
            realSizeReachS++;
        if (U.parentT > 0)
            realSizeReachT++;
    }
    assert(sizeS == realSizeS);
    assert(sizeT == realSizeT);
    assert(sizeReachS == realSizeReachS);
    assert(sizeReachT == realSizeReachT);
}


bool FlowCutter::checkIsCut(const std::vector<int>& cut, int ss, int tt) {
    std::vector<bool> isInCut(vertices.size(), false);
    std::vector<bool> visited(vertices.size(), false);
    for (const int v : cut)
        isInCut[v] = true;
    queue.clear();
    queue.push(ss);
    visited[s] = true;
    while (!queue.empty()) {
        int u = queue.pop();
        if (u == tt)
            return false;
        for (const int v : vertices[u].neigh) {
            if (visited[v])
                continue;
            if (isInCut[v])
                continue;
            visited[v] = true;
            queue.push(v);
        }
    }
    return true;
}