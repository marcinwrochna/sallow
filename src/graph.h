#pragma once
#include <cassert>
#include <functional>
#include <vector>

#include "utils.h"

typedef int VertexID;
typedef std::vector<VertexID> Ordering;


struct Bounds
{
    int lower, upper;
    void updateLower(int newLower) { assert(newLower <= upper); lower = std::max(lower, newLower); }
    void updateUpper(int newUpper) { assert(newUpper >= lower); upper = std::min(upper, newUpper); }
};


class BasicOneBasedGraph
{
 public:
    class VertexIDRange
    {
        const BasicOneBasedGraph& _graph;
     public:
        class const_iterator
        {
         public:
            constexpr explicit const_iterator(VertexID id) noexcept : _id(id) {}
            const_iterator(const const_iterator&) = default;
            const_iterator(const_iterator&&) = default;
            const_iterator& operator=(const const_iterator&) = default;
            const_iterator& operator=(const_iterator&&) = default;
            constexpr VertexID operator*() const noexcept { return _id; }
            constexpr const_iterator& operator++() noexcept { ++_id; return *this; }
            constexpr void operator++(int) noexcept { ++_id; }
            constexpr bool operator!=(const const_iterator& other) const noexcept { return _id != other._id; }
         private:
            VertexID _id;
        };
        constexpr explicit VertexIDRange(const BasicOneBasedGraph& graph) noexcept : _graph(graph) {}
        constexpr const_iterator begin() const noexcept { return const_iterator{1}; }
        constexpr const_iterator end() const noexcept { return const_iterator{1 + static_cast<VertexID>(_graph._nVertices)}; }
    };
    class NeighRange
    {
     public:
        typedef std::vector<VertexID>::iterator iterator;
        constexpr explicit NeighRange(std::vector<VertexID>& neigh) noexcept : _neigh(neigh) {}
        iterator begin() const noexcept { return _neigh.begin(); }
        iterator end() const noexcept { return _neigh.end(); }
        size_t size() const noexcept { return _neigh.size(); }
     private:
        std::vector<VertexID>& _neigh;
    };
    class ConstNeighRange
    {
     public:
        typedef std::vector<VertexID>::const_iterator const_iterator;
        constexpr explicit ConstNeighRange(const std::vector<VertexID>&  neigh) : _neigh(neigh) {}
        const_iterator begin() const noexcept { return _neigh.cbegin(); }
        const_iterator end() const noexcept { return _neigh.cend(); }
        size_t size() const noexcept { return _neigh.size(); }
     private:
        const std::vector<VertexID>&  _neigh;
    };

    BasicOneBasedGraph() : _nVertices(0), _adjList(1) {}
    explicit BasicOneBasedGraph(int nVertices) : _nVertices(nVertices), _adjList(nVertices + 1) {}
    BasicOneBasedGraph(BasicOneBasedGraph&& other) :
        _nVertices(other._nVertices),
        _adjList(std::move(other._adjList)) {}
    explicit BasicOneBasedGraph(std::vector<std::vector<VertexID>>&& adjList) :
        _nVertices(adjList.size() - 1),
        _adjList(std::move(adjList)) {
        assert(_adjList.size() >= 2);
        assert(_adjList[0].empty());
        for (int i = 1; i < _adjList.size(); ++i)
            assert(!contains(_adjList[i], 0));
    }

    constexpr VertexIDRange vertices() const noexcept { return VertexIDRange{*this}; }
    constexpr size_t nVertices() const noexcept { return _nVertices; }
    constexpr size_t maxVertexID() const noexcept { return _nVertices; }
    NeighRange operator[](VertexID id) noexcept { return NeighRange{_adjList[id]}; }
    ConstNeighRange operator[](VertexID id) const noexcept { return ConstNeighRange{_adjList[id]}; }
    void addEdge(VertexID u, VertexID v) noexcept { assert(!contains(_adjList[u], v)); _adjList[u].push_back(v); _adjList[v].push_back(u); }
    VertexID addVertex() noexcept { _adjList.emplace_back(); ++_nVertices; return _nVertices;  }
    void clear() noexcept { _adjList.clear(); _adjList.emplace_back(); _nVertices = 0;}
    void reset(int nVertices) { _adjList.clear(); _adjList.resize(nVertices + 1); _nVertices = nVertices;}

    static BasicOneBasedGraph getInducedSubgraph(const BasicOneBasedGraph& graph, const std::vector<VertexID> vertices) {
        std::vector<int> newID(graph.nVertices() + 1, 0);
        for (VertexID u = 1; u <= vertices.size(); ++u)
            newID[vertices[u - 1]] = u;
        std::vector<std::vector<VertexID>> result(vertices.size() + 1);
        for (VertexID u = 1; u <= vertices.size(); ++u) {
            for (int v : graph[vertices[u - 1]])
                if (newID[v])
                    result[u].push_back(newID[v]);
            std::sort(result[u].begin(), result[u].end());
        }
        return BasicOneBasedGraph{std::move(result)};
    }

 private:
    size_t _nVertices;
    std::vector<std::vector<VertexID>> _adjList;
};


template <typename OneBasedGraph, typename Predicate>
class InducedSubgraphView {
 public:
    InducedSubgraphView(const OneBasedGraph& graph, size_t size, Predicate predicate)
            : _graph(graph), _size(size), _predicate(predicate) {}
    constexpr size_t nVertices() const noexcept { return _size; }
    constexpr size_t maxVertexID() const noexcept { return _graph.maxVertexID(); }
    constexpr auto vertices() const noexcept { return FilterView{_graph.vertices(), _predicate}; }
    constexpr auto operator[](VertexID id) const noexcept { return FilterView{_graph[id], _predicate}; }
 private:
    const OneBasedGraph& _graph;
    const size_t _size;
    const Predicate _predicate;
};


template <typename OneBasedGraph>
class CutSubgraphView : public InducedSubgraphView<OneBasedGraph, std::function<bool(VertexID)>> {
 public:
    CutSubgraphView(const OneBasedGraph& graph, const std::vector<VertexID>& cut)
            : InducedSubgraphView<OneBasedGraph, std::function<bool(VertexID)>>(
                  graph,
                  graph.nVertices() - cut.size(),
                  [this](VertexID v) -> bool { return _inSubgraph[v]; }
              ),
              _inSubgraph(graph.maxVertexID() + 1, true) {
        for (VertexID v : cut)
            _inSubgraph[v] = false;
    }
 private:
    std::vector<bool> _inSubgraph;
};
