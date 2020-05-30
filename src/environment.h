#pragma once
#include <csignal>
#include <chrono>
#include <functional>
#include <random>

#include "tinyGraphDB.h"

/**
 * Global things needed by the algorithm (random generator, signal reference, etc.).
 */
class Environment
{
 public:
    typedef std::minstd_rand RandGen;
    typedef std::chrono::high_resolution_clock Clock;

    RandGen randGen;
    TinyGraphDB db;
    const volatile std::sig_atomic_t& signalled; // Set to true when the process is interrupted by a signal.
    std::chrono::time_point<Clock> constructedTime;
    std::chrono::time_point<Clock> startTime;

    Environment (unsigned int seed, const volatile std::sig_atomic_t& signalled_) :
        randGen(seed), signalled(signalled_),
        constructedTime(Clock::now()),
        startTime(constructedTime) {
        fprintf(stderr, "Initializing small graph database...\n");
        db.initialize();
    }

    void resetTime() { startTime = Clock::now(); }
    double msSinceReset() { return std::chrono::duration<double, std::milli>(Clock::now() - startTime).count(); }


    std::string benchmark(std::function<int()> f, int maxSeconds = 10, int minRepetitions = 3, int maxRepetitions = 1000) {
        resetTime();
        int reps = 0;
        long totalResult = 0;
        int minResult = 2'000'000'000;

        while (true) {
            int result = f();
            if (!result)
                result = 1'000'000'000;
            totalResult += result;
            minResult = std::min(minResult, result);
            reps++;
            if (signalled)
                break;
            if (reps >= maxRepetitions)
                break;
            if (reps >= minRepetitions && msSinceReset() > maxSeconds * 1000)
                break;
        }
        if (!reps)
            return "?";
        char buf[1000];
        std::snprintf(buf, sizeof(buf), "%ld (min %d) time %fms\t\t\treps %d", totalResult / reps, minResult, msSinceReset() / reps, reps);
        return std::string(buf);
    }
};
