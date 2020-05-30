#include "io.h"

#include <cassert>


BasicOneBasedGraph readGraph(std::istream& input) {
    BasicOneBasedGraph result;
    int nVertices, nEdges;
    int nReadEdges = 0;
    std::string line;
    assert(input);
	while (input) {
        std::getline(input, line);
        assert(!input.bad());
        // Skip empty and comment lines.
        if (line.empty() || line[0] == 'c')
            continue;
        // First non-comment line.
        else if (line[0] == 'p') {
            assert(result.nVertices() == 0);
            char dummy;
            std::string descriptor;
            std::istringstream iss{line};
            iss >> dummy >> descriptor >> nVertices >> nEdges;
            assert(!iss.fail());
            // assert(descriptor == "tdp");
            assert(nVertices >= 0);
            result.reset(nVertices);
        // Remaining lines: edges.
        } else {
            std::istringstream iss{line};
            VertexID u, v;
            iss >> u >> v;
            assert (1 <= u && 1 <= v);
            assert(u <= nVertices && v <= nVertices);
            if (u != v)
                result.addEdge(u, v);
            nReadEdges++;
            assert(!iss.fail());
        }
    }
    assert(!input.bad());
    assert(nEdges == nReadEdges);
    return result;
}


BasicOneBasedGraph readNautyGraph(std::istream& input) {
    BasicOneBasedGraph result;
    int nVertices, nEdges;
    int nReadEdges = 0;
    bool firstLine = true;
    std::string line;
    assert(input);
	while (input) {
        std::getline(input, line);
        assert(!input.bad());
        // Skip empty and header lines.
        if (line.empty() || line[0] == 'G')
            continue;
        std::istringstream iss{line};
        if (firstLine) {
            firstLine = false;
            iss >> nVertices >> nEdges;
            result.reset(nVertices);
        } else {
            //for (int i = 0; i < nEdges; ++i ) {
            while (iss.good()) {
                VertexID u, v;
                iss >> u >> v >> std::ws;
                assert (1 <= u && 1 <= v);
                assert(u <= nVertices && v <= nVertices);
                result.addEdge(u, v);
                nReadEdges++;
                assert(!iss.fail());
            }
        }
    }
    assert(nEdges == nReadEdges);
    return result;
}
