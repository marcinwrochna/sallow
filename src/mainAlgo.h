#pragma once
#include <memory>
#include <vector>

#include "environment.h"
#include "flowCutter.h"
#include "graph.h"
#include "graphUtils.h"
#include "simpleCutter.h"

class MainAlgo
{
    public:
        /** Graph should be connected. */
        MainAlgo(const BasicOneBasedGraph& graph_,
                 Environment& environment_,
                 int goodCutoff_ = 1,
                 int badCutoff_ = std::numeric_limits<int>::max(),
                 int recursionDepth_ = 0)
            : graph(graph_), environment(environment_),
              goodCutoff(goodCutoff_), badCutoff(badCutoff_),
              recursionDepth(recursionDepth_),
              doneStrength(-1), flowCutterCutoffs(0) {}

        /** Run the algorithm. */
        void run(int strength = 1'000'000'000);

        TreedepthResult bestResult;
    private:
        TreedepthResult recurseIntoComponents(const std::vector<VertexID>& cut, int strength);

        /** Find a cut using simpleCutter (smallest of minRetries) and recurse with cStrength. */
        bool doSimpleCut(int cStrength, int minRetries);
        /**
         * Find a cut of at least given balance using flowCutter and recurse.
         *
         * We start flowCutter with two random vertices;
         * they are replaced with most distant vertices `push` many times.
         * Thus push >= 2 will give similar cuts, push 0 or 1 gives more randomness.
         */
        bool doCut(double balance, int push, int cStrength);
        /** Update lower and upper bounds with a depth-first-search from s. */
        bool doDFS(VertexID s = 0);
        /** Update lower bound with degeneracy(graph) + 1. Only makes sense to call it once. */
        bool doDegeneracy();
        /** Run the superFastGreedy algorithm. */
        bool doSuperFastGreedy(int cutoff, int lookahead);
        /** Run the greedy elimination algorithm with coefficients alpha, beta, gamma. */
        bool doGreedy(int cutoff, std::tuple<int, int, int> coefficients);

        bool update(TreedepthResult&& other, const char* loggedName);
        bool updateLowerBound(int other, const char* loggedName);

        const BasicOneBasedGraph& graph;
        Environment& environment;
        int goodCutoff; // Stop as soon as we get treedepth at most this.
        int badCutoff; // Skip attempts that won't lead to treedepth strictly lower than this.
        int recursionDepth;
        int doneStrength; // Strength we've completed at the moment in `run()`.

        std::unique_ptr<FlowCutter> flowCutter;
        std::unique_ptr<SimpleCutter> simpleCutter;
        int flowCutterCutoffs; // How many times we did not find a small enough balanced cut.
};
