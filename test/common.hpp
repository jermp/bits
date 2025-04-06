#pragma once

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "external/doctest/doctest/doctest.h"

#include <iostream>
#include <vector>
#include <random>

#include "essentials.hpp"

namespace bits::test {

/*
    Generate [sequence_length] random integers, distributed as a poisson RV with a random mean.
*/
std::vector<uint64_t> get_sequence(const uint64_t sequence_length,                       //
                                   const uint64_t max_int = 1000,                        //
                                   const uint64_t seed = essentials::get_random_seed())  //
{
    srand(seed);
    double mean = rand() % (max_int + 1);
    std::mt19937 rng(seed);
    std::poisson_distribution<uint64_t> distr(mean);
    std::vector<uint64_t> seq(sequence_length);
    std::generate(seq.begin(), seq.end(), [&]() { return distr(rng) % (max_int + 1); });
    assert(seq.size() == sequence_length);
    return seq;
}

/*
    Generate a sorted sequence of [sequence_length] random integers,
    distributed as a poisson RV with a random mean.
*/
std::vector<uint64_t> get_sorted_sequence(const uint64_t sequence_length,                       //
                                          const uint64_t max_int = 1000,                        //
                                          const bool all_distinct = false,                      //
                                          const uint64_t seed = essentials::get_random_seed())  //
{
    srand(seed);
    double mean = rand() % (max_int + 1);
    std::mt19937 rng(seed);
    std::poisson_distribution<uint64_t> distr(mean);
    uint64_t universe = 0;
    std::vector<uint64_t> seq(sequence_length);
    std::generate(seq.begin(), seq.end(), [&]() {
        uint64_t val = distr(rng) % (max_int + 1);
        universe += val + all_distinct;
        return universe;
    });
    assert(seq.size() == sequence_length);
    assert(seq.back() == universe);
    assert(std::is_sorted(seq.begin(), seq.end()));
    return seq;
}

uint64_t get_random_uint(const uint64_t max_int = 1000,
                         const uint64_t seed = essentials::get_random_seed())  //
{
    srand(seed);
    return rand() % max_int + 1;  // at least 1
}

template <typename T>
void print(std::vector<T> const& v) {
    std::cout << '[';
    for (uint64_t i = 0; i != v.size(); ++i) {
        auto const& x = v[i];
        std::cout << x;
        if (i != v.size() - 1) std::cout << ',';
    }
    std::cout << ']';
    std::cout << std::endl;
}

}  // namespace bits::test

using namespace bits;
