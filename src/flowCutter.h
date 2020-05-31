#pragma once
#include <functional>
#include <vector>

#include "environment.h"
#include "graph.h"
#include "queue.h"


class FlowCutter
{
 public:
    struct Vertex {
        std::vector<int> neigh; // neighborhood list of this vertex, we copy it here to have faster access.
        bool inS, inT; // whether the vertex is in S/T (pierceS/T always is).
        // Invariant: S and T are disjoint and not connected by an edge.
        bool boundaryS, boundaryT; // whether the vertex is in S/T or in flowEndpointsS/T.
        // Invariant: vertices in flowEndpointsS/T are never in S/T.
        // Invariant: all neighbors of (S-pierceS) are boundaryS, some neighbors of pierceS are too.
        // (indeed we always have a maximum flow from (S-pierceS) to (T-pierceT), plus some flow from pierceS to pierceT.

        int inflow, outflow;  // predecessor/successor in flow; non-zero iff on some flow path.
        // Invariant: inflow[u] == u iff u is in flowEndpointsS. Same for outflow, T.
        // Invariant: u == inflow[v] iff outflow[u] == v for u != v.
        // Invariant: no vertices with inflow (equivalently, outflow) are in S nor T.

        int parentS, parentT;  // parent in BFS of residual digraph; non-zero iff reachable from from pierceS/T;
        // negative if only pre-reachable, see below. Note boundaryS has usually no parentS.
        // Invariant: inS/inT is true iff parentS/parentT points to itself.

       // The residual digraph from S has two copies pre-v and post-v of every original vertex v, and:
        // - an arc from pre-v to post-v for all v not in flow, vice-versa for all v in flow;
        // - an arc from post-u to pre-v for all arcs u->v not in flow, vice-versa for all u->v in flow.
        // (Here if u->v is in flow then v->u is not, so the new arcs go back the flow only).
        // If a vertex is post-reachable (i.e. pre-v is reachable from S) then it is also pre-reachable,
        // because either its not in flow, in which case its only accessible via pre-v,
        // or its in the flow, in which case there is an arc from post-v to pre-v.
        // Note that an augmenting path can go through pre-v and then much later through post-v.
        // Note also that v may be pre-reachable from S and from T (i.e. parentS,parentT < 0)
        // without giving an augmenting path.
    };
    enum class EndReason { none, complete, cutoff, signal, badST, callback };
    typedef std::function<EndReason(std::vector<VertexID>, EndReason, int)> YieldCutCallback;

 private:
    const BasicOneBasedGraph& graph;
    Environment& environment;
    std::vector<Vertex> vertices; // Map from vertex ID to Vertex. Index 0 is dummy, vertex IDs are 1-based.
    std::vector<int> distS, distT; // Original distance from first s/t vertex. Rarely needed, so we store it outside of vertices.

    int s, t;
    int pierceS, pierceT;
    int sizeS, sizeT; // Number of nodes in S/T.
    int sizeReachS, sizeReachT;  // Number of nodes with positive parentS/T. Only valid after bfsGrowParentS/T returns 0.
    std::vector<int> flowEndpointsS;
    std::vector<int> flowEndpointsT;
    bool clearedParentS, clearedParentT; // Whether parentS, parentT are cleared (i.e. parentS=0 if not in S).
    Queue<int> queue; // Queue reused by all BFSs.
public:
    YieldCutCallback yieldCutCallback; // Called in yieldCutS/yieldCutT.

    int cutoff; // We already have a treedepth decomposition of this depth, stop if we can't improve it,
    // e.g. when cut size would increase to cutoff-1.

    /** Initialize with a graph with vertices 1..n (0 is dummy). */
    FlowCutter(const BasicOneBasedGraph&, Environment&);

    EndReason findCuts(int _s, int _t, int _cutoff, YieldCutCallback);

    /** Update heuristics and return whether the next pierce side should be S, rather than T. */
    bool choosePierceSide();

    /** Called whenever we find a new cut. There might be more balanced cuts of equal size later. */
    void yieldSoftCut(std::vector<int> cut);
    /** Called whenever we find a new cut and all future cuts will be larger. Return whether to stop. */
    EndReason yieldCutS(EndReason r);
    EndReason yieldCutT(EndReason r);

 private:
    void addToS(int v);
    void addToT(int v);
    void pushFlowEndpointS(int index);
    void pushFlowEndpointT(int index);

    /** Grow BFS from pierceS. Return boundary of T if found (giving an ST-path), zero otherwise. */
    int bfsGrowParentS();
    /** Grow BFS from pierceT. Return boundary of S if found (giving an ST-path), zero otherwise. */
    int bfsGrowParentT();

    void bfsGrowInS();
    void bfsGrowInT();

    /** Augment the flow given a neighbor of pierceT returned by bfsGrowParentS. Clears parentS/T. */
    void augmentFlowViaS(int last);
    /** Augment the flow given a neighbor of pierceS returned by bfsGrowParentT. Clears parentS/T. */
    void augmentFlowViaT(int first);

    void clearParentS();
    void clearParentT();

    /*  Selects a new pierceS among flowEndpointsS (assuming bfsGrowInS was just called).
        Returns true if we want to stop.
    */
    EndReason selectNewPierceS();
    EndReason selectNewPierceT();

    /** Initialize distS/T by BFS from s/t. */
    void initDistS();
    void initDistT();

    void checkInvariants();
    bool checkIsCut(const std::vector<int>& cut, int ss, int tt);

 public:
    class ReachableSubgraphViewS {
    private:
        const BasicOneBasedGraph& _graph;
        const FlowCutter& _flowCutter;
    public:
        ReachableSubgraphViewS(const BasicOneBasedGraph& graph, const FlowCutter& flowCutter)
                : _graph(graph), _flowCutter(flowCutter) {}
        constexpr auto vertices() const noexcept {
            return FilterView{_graph.vertices(), [this](VertexID v) -> bool { return _flowCutter.vertices[v].parentS > 0; } };
        }
        constexpr size_t nVertices() const noexcept { return _flowCutter.sizeReachS; }
        constexpr size_t maxVertexID() const noexcept { return _graph.maxVertexID(); }
        auto operator[](VertexID id) const noexcept {
            return FilterView{_graph[id], [this](VertexID v) -> bool { return _flowCutter.vertices[v].parentS > 0; } };
        }
    };
   class ReachableSubgraphViewT {
    private:
        const BasicOneBasedGraph& _graph;
        const FlowCutter& _flowCutter;
    public:
        ReachableSubgraphViewT(const BasicOneBasedGraph& graph, const FlowCutter& flowCutter)
                : _graph(graph), _flowCutter(flowCutter) {}
        constexpr auto vertices() const noexcept {
            return FilterView{_graph.vertices(), [this](VertexID v) -> bool { return _flowCutter.vertices[v].parentT > 0; } };
        }
        constexpr size_t nVertices() const noexcept { return _flowCutter.sizeReachT; }
        constexpr size_t maxVertexID() const noexcept { return _graph.maxVertexID(); }
        auto operator[](VertexID id) const noexcept {
            return FilterView{_graph[id], [this](VertexID v) -> bool { return _flowCutter.vertices[v].parentT > 0; } };
        }
    };
};
