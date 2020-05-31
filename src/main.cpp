#include <cassert>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <chrono>

#define RELEASE

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
    #ifndef RELEASE
        filename = "data/pace20-td-heur-pub/heur_001.gr";
    #endif
    if (argc > 1 && std::strlen(argv[1]) > 0)
        filename = argv[1];

    int goodCutoff = 0;
    #ifndef RELEASE
        // Best results known for public tests; to speed up local testing we stop as soon as we get it.
        const std::array<int, 100> optil{45,14,38,12,34,11,19,60,18,7,19,138,15,17,10,17,28,25,54,8,45,18,28,175,17,17,47,14,60,15,26,59,533,106,138,45,96,30,20,448,31,744,37,11,132,154,891,11,219,1262,44,146,770,901,1338,2286,1434,45,12,381,332,87,90,12,69,103,768,13,58,5694,669,14,63,167,137,3072,162,173,68,163,16,194,158,253,214,148,212,162,3,293,34439,274,10990,18,390,17850,1745,77758,394,397};
        if (filename.length() >= 11 && filename.substr(filename.length() - 11, 5) == "heur_")
            goodCutoff = optil[std::stoi(filename.substr(filename.length() - 6, 3)) / 2];
        fprintf(stderr, "goodCutoff = %d\n", goodCutoff);
    #endif

    fprintf(stderr, "%s:\t", filename.c_str());
    std::istream* input = &std::cin;
    std::ifstream ifstream;
    if (filename != "-") {
        ifstream.open(filename, std::ios::binary);
        input = &ifstream;
    }
    BasicOneBasedGraph graph = readGraph(*input);
    fprintf(stderr, "%zu vertices\n", graph.nVertices());
    // int n = 11;
    // BasicOneBasedGraph graph(n);
    // graph.addEdge(1, 2);
    // graph.addEdge(2, 3);
    // graph.addEdge(3, 4);
    // graph.addEdge(4, 5);
    // graph.addEdge(1, 6);
    // graph.addEdge(6, 7);
    // graph.addEdge(7, 3);
    // graph.addEdge(3, 10);
    // graph.addEdge(10, 11);
    // graph.addEdge(11, 5);
    // graph.addEdge(2, 8);
    // graph.addEdge(8, 9);
    // graph.addEdge(9, 4);

    MainAlgo mainAlgo{graph, environment, goodCutoff};
    mainAlgo.run();
    #ifdef RELEASE
        if (mainAlgo.bestResult.depth)
            printDecomposition(mainAlgo.bestResult);
    #else
        printf("%d\n", mainAlgo.bestResult.depth);
        fprintf(stderr, "goodCutoff = %d\n", goodCutoff);
    #endif

    return EXIT_SUCCESS;
}
