
#pragma once
#include <sstream>

#include "graph.h"

/** Read a graph in DIMACS format. */
BasicOneBasedGraph readGraph(std::istream& input);

/** Read a single graph as outputted by `nauty-showg -e -p1 -o1`. */
BasicOneBasedGraph readNautyGraph(std::istream& input);
