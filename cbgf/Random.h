//
// Created by Xun Zheng on 10/10/16.
//

#ifndef CBGF_THREADLOCALRANDOM_H
#define CBGF_THREADLOCALRANDOM_H


#include <random>
#include <thread>
#include <chrono>

/**
 * Random utilities.
 */
class Random {
public:
    /**
     * Construct random engine from thread-safe seeds.
     */
    Random() {
        std::hash<std::thread::id> hasher;
        auto clock = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        rng_.seed((unsigned int) (hasher(std::this_thread::get_id()) + clock));
    }

    /**
     * Draw an element uniformly from {0, 1, ..., k-1}.
     *
     * @param k exclusive upper bound
     * @return a random integer from [0, k).
     */
    int next_int(int k) {
        return (int) (uniform_(rng_) * k);
    }

    /**
     * Toss a fair coin.
     *
     * @return a random boolean
     */
    bool next_bool() {
        return bernoulli_(rng_);
    }

    /**
     * Draw a real number uniformly from [0, 1).
     *
     * @return a random double from uniform(0,1)
     */
    double next_double() {
        return uniform_(rng_);
    }

private:
    std::mt19937 rng_;
    std::bernoulli_distribution bernoulli_;
    std::uniform_real_distribution<double> uniform_;
};


#endif //CBGF_THREADLOCALRANDOM_H
