#include <cassert>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <chrono>

#include "environment.h"
#include "utils.h"
#include "io.h"
#include "graph.h"
#include "graphUtils.h"
#include "mainAlgo.h"


volatile std::sig_atomic_t signalled;
extern "C" void signal_handler(int) { signalled = 1; }


static void printDecomposition(const TreedepthResult& result) {
    printf("%d\n", result.depth);
    for (VertexID u = 1; u < result.parent.size(); ++u)
        printf("%d\n", result.parent[u]);
}


int main(int argc, char* argv[]) {
    if (argc > 1 && argv[1] == std::string("--version")) {
        printf("v4 submitted\n");
        return EXIT_SUCCESS;
    }

    signalled = 0;
    std::signal(SIGTERM, signal_handler);
    std::signal(SIGINT, signal_handler);
    std::signal(SIGPIPE, signal_handler);

    unsigned int seed = std::random_device()();
    Environment environment(seed, signalled);
    fprintf(stderr, "random seed = %u\n", seed);

    std::string filename = "-";
    //     "data/pace20-td-heur-pub/heur_001.gr";
    if (argc > 1 && std::strlen(argv[1]) > 0)
        filename = argv[1];
    fprintf(stderr, "%s:\t", filename.c_str());
    std::istream* input = &std::cin;
    std::ifstream ifstream;
    if (filename != "-") {
        ifstream.open(filename, std::ios::binary);
        input = &ifstream;
    }

    BasicOneBasedGraph graph = readGraph(*input);
    fprintf(stderr, "%zu vertices\n", graph.nVertices());

    MainAlgo mainAlgo{graph, environment};
    mainAlgo.run();
    if (mainAlgo.bestResult.depth)
        printDecomposition(mainAlgo.bestResult);

    return EXIT_SUCCESS;
}
