#include "mainAlgo.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <string>

#include "greedy.h"
#include "utils.h"


namespace {
    // Coefficients to use in doGreedy if we execute it before/after attempting doCuts.
    constexpr const std::array<std::tuple<int, int, int>, 14> preCoeffs = {{
        {1, 6, 0}, {1, 2, 0}, {4, 3, 0}, {0, 9, 1},
        {1, 4, 0}, {2, 3, 0}, {0, 1, 0}, {1, 1, 0},
        {1, 20, 0}, {3, 8, 0}, {5, 50, 1}, {0, 40, 1},
        {1, 3, 0}, {1, 5, 0}
    }};
    constexpr const std::array<std::tuple<int, int, int>, 3> postCoeffs = {{
        {3, 50, 1}, {1, 50, 5}, {1, 100, 1}
    }};
}


void MainAlgo::run(int strength) {
    int nVertices = graph.nVertices();
    assert(nVertices > 0);
    if (nVertices <= 10)
        strength = std::min(strength, 6);

    // Strength 0.
    if (doDFS() || environment.signalled)
        return;
    doneStrength = 0;
    if (doneStrength == strength)
        return;

    // Strength 1.
    int initialCutoff = 5 * badCutoff + 5;
    if (doDFS(1) || environment.signalled)
        return;
    if (doSuperFastGreedy(std::min(bestResult.depth, initialCutoff + 1), 2) || environment.signalled)
        return;
    doneStrength = 1;
    if (doneStrength == strength)
        return;

    // Strength 2.
    if (doDegeneracy() || environment.signalled)
        return; // This takes care of cliques in particular.
    if (doSuperFastGreedy(std::min(bestResult.depth, initialCutoff + 1), 64) || environment.signalled)
        return;
    for (int i = 0; i < bitlength(nVertices); ++i)
        if (doDFS() || environment.signalled)
            return;
    doneStrength = 2;
    if (doneStrength == strength)
        return;
    if (bestResult.depth > initialCutoff && nVertices < 5000) // If our best attempt at this point is so bad, give up.
        return;

    // Strength 3.
    for (int i = 0; i < bitlength(nVertices); ++i)
        if (doDFS() || environment.signalled)
            return;
    // How many greedy runs to perform at first (strength 3 & 4): avoid too large, too deep graphs.
    int nGreedies = preCoeffs.size();
    if (bestResult.depth > 100)
        nGreedies = std::min(nGreedies, strength + 3);
    if ((nVertices * (long) bestResult.depth / 1'000'000) > 20'000)
        nGreedies = 0;
    if (nVertices > 150'000)
        nGreedies = std::min(nGreedies, recursionDepth ? 1 : 2);
    if ((nVertices * (long) bestResult.depth / 1'000'000) > 40 && recursionDepth)
        nGreedies = std::min(nGreedies, 1);
    log(recursionDepth, "nGreedies = %d factor %ld\n", nGreedies, (nVertices * (long) bestResult.depth / 1'000'000));
    // We don't expect to get results worse than OPT / expectedRatio.
    double expectedRatio;
    if (nVertices < 10'000)
        expectedRatio = 0.5;
    else
        expectedRatio = 0.4;

    if (nGreedies)
        if (doGreedy(std::min(bestResult.depth, (int)(badCutoff / expectedRatio) + 2), {1, 9, 0}) || environment.signalled)
            return;
    if (doCut(0.4, 3, 3) || environment.signalled)
            return;
    doneStrength = 3;
    if (doneStrength == strength)
        return;
    if (bestResult.depth > badCutoff / expectedRatio + 1)
        return;

    // Strength 4.
    if (nVertices < 100)
        expectedRatio = 0.75;
    else if (nVertices < 250'000)
        expectedRatio = 0.65;
    else
        expectedRatio = 0.5;
    for (int i = 0; i < nGreedies - 1; ++i)
        if (doGreedy(std::min(bestResult.depth, (int)(badCutoff / expectedRatio) + 2), preCoeffs[i]) || environment.signalled)
            return;
    double balance = std::uniform_real_distribution<>(0.34, 0.5)(environment.randGen);
    int push = std::uniform_int_distribution<>(0,1)(environment.randGen);
    if (doCut(balance, 1, 4) || environment.signalled)
            return;
    doneStrength = 4;
    if (doneStrength == strength)
        return;
    if (bestResult.depth > badCutoff / expectedRatio + 1)
        return;

    // Strength 5.
    for (int i = 0; i < std::min(nGreedies - 2, (int)postCoeffs.size()); ++i)
          if (doGreedy(bestResult.depth, postCoeffs[i]) || environment.signalled)
              return;
    for (int rep = 0; rep <= strength - 5; ++rep) {
        std::uniform_real_distribution<> balanceDistribution(0.3, 0.54);
        std::uniform_int_distribution<> pushDistribution(0, 1);
        balance = balanceDistribution(environment.randGen);
        push = pushDistribution(environment.randGen);
        int substrength = 4;
        if (recursionDepth)
            substrength = rep / 2 + 4;
        else
            substrength = std::min(std::max(4, rep / 3 + 4), 6);
        if (rep % 10 == 9) {
            if (doSimpleCut(substrength, 1) || environment.signalled)
                return;
        } else {
            if (doCut(balance, push, substrength) || environment.signalled)
                return;
        }
        if (nGreedies > 2 || (nGreedies > 1 && rep % 6 == 5))
            if (doGreedy(badCutoff, (rep % 3) ? preCoeffs[rep % 7] : postCoeffs[rep % postCoeffs.size()]) || environment.signalled)
                return;
        doneStrength = rep + 5;
    }
    /*{
        bestResult.computeHeights();
        bestResult.computeDepths();
        auto s = benchmark([&](){ return  dfsPreOrder(graph, getRandomVertex(graph, environment.randGen)).depth; }, 1);
        log(recursionDepth - 1, "\t\t\t\t\t\tDFS = %s\n", s.c_str());
        s = benchmark([&](){ return  degeneracy(graph); }, 1);
        log(recursionDepth - 1, "\t\t\t\t\t\tdegeneracy = %s\n", s.c_str());
        for (int lookahead : {2, 64}) {
            s = benchmark([&](){ return superFastGreedy(graph, bestResult.ordering, lookahead, nVertices + 1, environment).depth; }, 1);
            log(recursionDepth - 1, "\t\t\t\t\t\tsuperFastGreedy(%d) = %s\n", lookahead, s.c_str());
        }
        s = benchmark([&](){ return greedyElimByHeur(graph, bestResult.heights, nVertices + 1, 1, 9, 0, environment).depth; });
        log(recursionDepth - 1, "\t\t\t\t\t\tgreedyElim(%d %d %d) = %s\n", 1, 9, 0, s.c_str());
    }*/
}


TreedepthResult MainAlgo::recurseIntoComponents(const std::vector<VertexID>& cut, int strength) {
    TreedepthResult result;
    if (cut.size() + 1 >= badCutoff) {
        result.clear();
        return result;
    }
    result.depth = 0; // Will be maximum component depth + cut.size().
    result.ordering.reserve(graph.nVertices());
    result.parent = std::vector<VertexID>(graph.maxVertexID() + 1, -1);
    for (VertexID v : cut)
        result.ordering.push_back(v);
    VertexID parent = 0;
    for (VertexID v : cut) {
        result.parent[v] = parent;
        parent = v;
    }
    std::vector<std::vector<VertexID>> bigComponents;
    std::vector<bool> inCut(graph.maxVertexID() + 1, false);
    for (VertexID v : cut)
        inCut[v] = true;
    std::vector<bool> visited = inCut;
    std::vector<VertexID> stack;
    for (VertexID s : graph.vertices()) {
        if (visited[s])
            continue;
        std::vector<VertexID> currentComponent;
        int nSubEdges = 0; // Twice the number of edges.
        stack.push_back(s);
        visited[s] = true;
        while (!stack.empty()) {
            VertexID u = stack.back();
            stack.pop_back();
            currentComponent.push_back(u);
            for (VertexID v : graph[u]) {
                if (!inCut[v])
                    ++nSubEdges;
                if (visited[v])
                    continue;
                stack.push_back(v);
                visited[v] = true;
            }
        }
        assert(nSubEdges % 2 == 0);
        int nSubVertices = currentComponent.size();
        assert(nSubVertices + 1 + cut.size() <= graph.nVertices());
        if (nSubEdges == nSubVertices * (nSubVertices - 1)) { // Immediately handle small cliques.
            result.depth = std::max(result.depth, nSubVertices);
            result.lowerBound = std::max(result.lowerBound, nSubVertices);
            for (VertexID v: currentComponent)
                result.ordering.push_back(v);
            VertexID subParent = parent;
            for (VertexID v: currentComponent) {
                result.parent[v] = subParent;
                subParent = v;
            }
        } else if (nSubVertices <= environment.db.maxSize) {
            const TinyTreedepthResult& r = environment.db.getDecompositionConnected(graph, currentComponent);
            result.depth = std::max(result.depth, (int)r.depth);
            result.lowerBound = std::max(result.lowerBound, (int)r.depth);
            for (int v : r.ordering)
                result.ordering.push_back(currentComponent[v]);
            for (int v = 0; v < currentComponent.size(); ++v)
                if (r.parent[v] != -1)
                    result.parent[currentComponent[v]] = currentComponent[r.parent[v]];
                else
                    result.parent[currentComponent[v]] = parent;
        } else {
            bigComponents.push_back(std::move(currentComponent));
        }
    }
    std::sort(bigComponents.begin(), bigComponents.end(), [](auto& a, auto& b){ return a.size() > b.size(); });
    for (auto currentComponent : bigComponents) {
        size_t nSubVertices = currentComponent.size();
        auto subgraph = BasicOneBasedGraph::getInducedSubgraph(graph, currentComponent);
        int subGoodCutoff = goodCutoff - (int)cut.size();
        subGoodCutoff = std::max(subGoodCutoff, result.depth); // Cut subalgorithm if it gets at least as good as max of other components.
        if (doneStrength >= 3) // Or if it gets better than what's needed to improve by 60%/40%/20%.
            subGoodCutoff = std::max(subGoodCutoff, (int)((badCutoff - (int)cut.size()) * 0.4 - 3));
        if (doneStrength >= 4)
            subGoodCutoff = std::max(subGoodCutoff, (int)((badCutoff - (int)cut.size()) * 0.6 - 3));
        if (doneStrength >= 5)
            subGoodCutoff = std::max(subGoodCutoff, (int)((badCutoff - (int)cut.size()) * 0.8 - 3));
        MainAlgo subAlgo{subgraph, environment, subGoodCutoff, badCutoff - (int)cut.size(), recursionDepth + 1};
        if (nSubVertices > 10)
            log(recursionDepth, "Begin(%d) %zu vertices\n", strength, nSubVertices);
        subAlgo.run(strength);
        const TreedepthResult& subResult = subAlgo.bestResult;
        if (nSubVertices > 10)
            log(recursionDepth, "End(%d) %zu vertices %d<=td<=%d (should be < %d)\n", strength, nSubVertices, subResult.lowerBound,
                subResult.depth, badCutoff - (int)cut.size());
        if (!subResult.depth || subResult.depth >= badCutoff - cut.size()) {
            result.clear();
            return result;
        }
        result.depth = std::max(result.depth, subResult.depth);
        result.lowerBound = std::max(result.lowerBound, subResult.lowerBound);
        for (VertexID v: subResult.ordering)
            result.ordering.push_back(currentComponent[v - 1]);
        for (int v = 1; v <= nSubVertices; ++v)
            if (subResult.parent[v])
                result.parent[currentComponent[v - 1]] = currentComponent[subResult.parent[v] - 1];
            else
                result.parent[currentComponent[v - 1]] = parent;
    }

    result.depth += cut.size();
    if (doneStrength >= 3) {
        TreedepthResult greedy = superFastGreedy(graph, result.ordering, 64, badCutoff, environment);
        if (greedy.depth && greedy.depth < result.depth) {
            log(recursionDepth, "Imp greedy on cut %d to %d\n", result.depth, greedy.depth);
            result = std::move(greedy);
        }
    }
    return result;
}


/**
 * balance - at what size of (side + cut) we should use the cut.
 * push - how many times we replace s/t with mostDistant(t/s).
 * cStrength - passed to runIteration on each component after cut.
 * bStrength - passed to runIteration whenever we update bounds in a callback from flowCutter.
 */
bool MainAlgo::doCut(double balance, int push, int cStrength) {
    if (flowCutterCutoffs >= 3 && bestResult.depth > 1500)
        return doSimpleCut(cStrength, 3);
    log(recursionDepth, "Running cut(balance %f, push %d, cStrength %d)...\n", balance, push, cStrength);
    if (!flowCutter)
        flowCutter = std::make_unique<FlowCutter>(graph, environment);
    std::vector<int> cut;
    for (int retry = 1; retry <= 20; ++retry) {
        int s = getRandomVertex(graph, environment.randGen);
        int t = getRandomVertex(graph, environment.randGen);
        for (int i = 0; i < push; ++i)
            if (i % 2)
                t = mostDistant(graph, s);
            else
                s = mostDistant(graph, t);
        if (s == t)
            continue;
        for (int v : graph[s])
            if (v == t)
                continue;

        int prevSize = 0;
        FlowCutter::YieldCutCallback yieldCallback =
            [this, &cut, balance, &prevSize](std::vector<VertexID> newCut, FlowCutter::EndReason r, int sizeST) {
                if (newCut.size() > 1.5 * prevSize) {
                    log(recursionDepth, "yield cut %zu + sizeST %d is %.3lf%%\tdiff=%d\n", newCut.size(), sizeST, 100*static_cast<double>(sizeST + newCut.size()) / graph.nVertices(), sizeST - (int)newCut.size());
                    prevSize = newCut.size();
                }
                if (sizeST + newCut.size() >= balance * graph.nVertices())
                    r = FlowCutter::EndReason::callback;
                if (r != FlowCutter::EndReason::none) {
                    cut = newCut;
                    log(recursionDepth, "stop cut %zu + sizeST %d is %.3lf%%\tdiff=%d\n", newCut.size(), sizeST, 100*static_cast<double>(sizeST + newCut.size()) / graph.nVertices(), sizeST - (int)newCut.size());
                }
                return r;
        };
        int cutoff = badCutoff;
        if (graph.nVertices() > 150'000)
            cutoff = std::min(cutoff, 3000);
        if (badCutoff >= 10'000)
            cutoff = std::min(cutoff, 300);
        FlowCutter::EndReason reason = flowCutter->findCuts(s, t, cutoff, yieldCallback);
        if (environment.signalled)
            return true;
        if (reason == FlowCutter::EndReason::cutoff) {
            ++flowCutterCutoffs;
            return doSimpleCut(cStrength, 3);
        } else if (reason == FlowCutter::EndReason::badST) {
            continue; // Retry.
        }
        if (!cut.empty() && cut.size() + 1 < badCutoff)
            break;
    }
    if (cut.empty())
        return doSimpleCut(cStrength, 3);
    log(recursionDepth, "Cut %zd\n", cut.size());
    return update(recurseIntoComponents(cut, cStrength), "cut");
}


bool MainAlgo::doSimpleCut(int cStrength, int minRetries) {
    log(recursionDepth, "Running simpleCut(cStrength %d, minRetries %d)...\n", cStrength, minRetries);
    if (!simpleCutter)
        simpleCutter = std::make_unique<SimpleCutter>(graph, environment);
    std::vector<VertexID> bestCut;
    int retry = 0;
    while ((bestCut.empty() || retry < minRetries) && retry < 500) {
        ++retry;
        if (environment.signalled)
            return true;
        std::vector<VertexID> cut = simpleCutter->getGrowingCut(badCutoff);
        if (!cut.empty() && cut.size() < badCutoff && (!bestCut.size() || cut.size() < bestCut.size()))
            bestCut = std::move(cut);
    }
    log(recursionDepth, "Simple cut %zd (bad %d best %d, retry %d)\n", bestCut.size(), badCutoff, bestResult.depth, retry);
    if (bestCut.empty())
        return false;
    return update(recurseIntoComponents(bestCut, cStrength), "simpleCut");
}


bool MainAlgo::doDFS(VertexID s) {
    if (!s)
        s = getRandomVertex(graph, environment.randGen);
    if (!recursionDepth && s == 1)
        log(recursionDepth, "Running DFS...\n");
    return update(dfsPreOrder(graph, s), "DFS");
}


bool MainAlgo::doDegeneracy() {
    log(recursionDepth, "Running degeneracy...\n");
    return updateLowerBound(degeneracy(graph) + 1, "degeneracy");
}


bool MainAlgo::doSuperFastGreedy(int cutoff, int lookahead) {
    log(recursionDepth, "Running superFastGreedy(%d) cutoff %d...\n", lookahead, cutoff);
    TreedepthResult r = superFastGreedy(graph, bestResult.ordering, lookahead, cutoff, environment);
    return update(std::move(r), "superFastGreedy");
}


bool MainAlgo::doGreedy(int cutoff, std::tuple<int, int, int> coefficients) {
    auto [alpha, beta, gamma] = coefficients;
    log(recursionDepth, "Running greedy(%d, %d, %d) cutoff %d...\n", alpha, beta, gamma, cutoff);
    bestResult.computeHeights();
    TreedepthResult r = greedyElimByHeur(graph, bestResult.heights, cutoff, alpha, beta, gamma, environment);
    return update(std::move(r), "greedy");
}


bool MainAlgo::update(TreedepthResult&& other, const char* loggedName) {
    other.assertCorrect(graph);
    int otherLowerBound = other.lowerBound;
    if (other.depth && (!bestResult.depth || other.depth < bestResult.depth)) {
        if (doneStrength >= 2) {
            while (true) {
                TreedepthResult local = locallyImprove(graph, other, environment);
                if (local.depth && local.depth < other.depth)
                    other = std::move(local);
                else
                    break;
            }
        }
        if (!recursionDepth)
            log(-1, "Improved %zd vertices: %d<=td<=%d (was %d) by %s.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
            graph.nVertices(), bestResult.lowerBound, other.depth, bestResult.depth, loggedName);
        else if (graph.nVertices() > 100)
            log(recursionDepth, "Improved %zd vertices: %d<=td<=%d (was %d) by %s.-------------\n",
            graph.nVertices(), bestResult.lowerBound, other.depth, bestResult.depth, loggedName);
        other.lowerBound = bestResult.lowerBound;
        bestResult = std::move(other);
        badCutoff = std::min(bestResult.depth, badCutoff);
    }
    if (bestResult.depth <= goodCutoff)
        log(recursionDepth - 1, "Goodcutoff\n");
    return updateLowerBound(otherLowerBound, loggedName) || bestResult.depth <= goodCutoff;
}


bool MainAlgo::updateLowerBound(int other, const char* loggedName) {
    if (other > bestResult.lowerBound)
        log(recursionDepth, "Lower bound improved %zd vertices: %d<=td<=%d (was %d) by %s.\n",
            graph.nVertices(), other, bestResult.depth, bestResult.lowerBound, loggedName);
    bestResult.lowerBound = std::max(bestResult.lowerBound, other);
    assert(bestResult.lowerBound <= bestResult.depth);
    assert(bestResult.depth >= badCutoff);
    return bestResult.lowerBound >= badCutoff;
}
