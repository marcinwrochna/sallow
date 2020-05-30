#pragma once
#include "environment.h"
#include "graph.h"
#include "graphUtils.h"


/** Return the degeneracy of a graph. */
int degeneracy(const BasicOneBasedGraph& graph);

/**
 * Greedy elimination by ordering, but checking the height of `lookahead` next vertices.
 *
 * Lookahead 2 is a very nice fast default and often the same order of magnitude as much larger.
 * Lookahead 1 is only rarely faster, only possibly worth it for n > 250'000.
 * Lookahead 64 is usually worth it for dense graphs, best depth >= 1000.
 * Larger lookaheads like 256 are usually either not significantly better than 64, or
 * slower than fastGreedyElim. Only possibly worth it when best depth > 100'000.
 *
 * Returns a null decomposition (depth 0) if we couldn't get depth < cutoff.
 */
TreedepthResult superFastGreedy(const BasicOneBasedGraph& graph, const Ordering& ordering,
                           int lookahead, int cutoff,
                           Environment& environment);

/**
 * Greedy elimination by alpha * degree (after contractions) + beta * height + gamma * centrality.
 *
 * Centrality may be any heuristic that should be larger for higher nodes in a treedepth decomposition.
 * A good default is to use alpha=1, beta=9, gamma=0, but trying diverse parameters helps.
 * Returns a null decomposition (depth 0) if we couldn't get depth < cutoff.
 */
TreedepthResult greedyElimByHeur(const BasicOneBasedGraph& graph,
                                 const std::vector<int>& centralities,
                                 int cutoff,
                                 int alpha, int beta, int gamma,
                                 Environment& environment);

/** Recomputes a tree decomposition bottom-up by trying to swap a node with its parent. */
TreedepthResult locallyImprove(const BasicOneBasedGraph& graph, TreedepthResult& old, Environment& environment);

/** Return the shortest prefix whose removal leaves at least two connected components of size >= minSize. */
std::vector<VertexID> cutFromOrdering(const BasicOneBasedGraph& graph, const Ordering& ordering, int minSize);

/**
 * Build a tree decomposition top down.
 *
 * Each vertex is added below its lowest processed neighbor, if the neighbors form a chain;
 * otherwise it is added above their least common ancestor.
 * This heuristic doesn't seem to ever work well.
 */
TreedepthResult topDownGreedy(const BasicOneBasedGraph& graph, const Ordering& ordering);

/**
 * Lower bound based on minimum degree of a minor.
 *
 * Contracts adjacent vertices of minimum degree to try getting as dense a minor as possible.
 * It's also a lower bound on treewidth, might be more useful for it.
 * For treedepth it's better than degeneracy, but still ~2x less than best upper bound, so hardly useful.
 */
int greedyMinorBound(const BasicOneBasedGraph& graph, int prevLowerBound, Environment& environment);
