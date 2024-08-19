#include <iostream>
#include <random>  // for poisson_distribution
#include <vector>

#include "include/elias_fano.hpp"

using namespace bits;

int main(int argc, char const** argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " [sequence_length]" << std::endl;
        return 1;
    }

    uint64_t n = std::stoull(argv[1]);  // sequence length
    uint32_t universe = 0;

    /* generate some n random integers, distributed as a poisson RV with a random mean. */
    uint64_t seed = essentials::get_random_seed();
    srand(seed);
    double mean = rand() % 1000;
    std::cout << "generating " << n << " integers destributed as Poisson(mean=" << mean << ")..."
              << std::endl;
    std::mt19937 rng(seed);
    std::poisson_distribution<uint32_t> distr(mean);

    std::vector<uint32_t> vec(n);
    std::generate(vec.begin(), vec.end(), [&]() {
        uint32_t val = distr(rng);
        universe += val;
        return universe;
    });

    /* another small test instance */
    // std::vector<uint32_t> vec{0, 3, 3, 5, 6, 6, 6, 9, 13, 23};
    // n = vec.size();
    // universe = vec.back();

    // std::vector<uint32_t> vec{0, 0, 0, 0, 0, 1, 1, 1, 1, 6, 6, 9, 13, 23, 23, 23, 45, 45, 88,
    // 88}; n = vec.size(); universe = vec.back();

    // std::vector<uint32_t> vec{3, 3, 6, 6, 9, 13, 23, 23, 23, 45, 45, 88, 88};
    // n = vec.size();
    // universe = vec.back();

    assert(vec.back() == universe);
    assert(std::is_sorted(vec.begin(), vec.end()));

    constexpr bool index_zeros = true;
    constexpr bool encode_prefix_sum = false;

    std::cout << "encoding them with elias_fano<index_zeros=" << index_zeros
              << ", encode_prefix_sum=" << encode_prefix_sum << ">..." << std::endl;
    elias_fano<index_zeros, encode_prefix_sum> ef;
    ef.encode(vec.begin(), vec.size());

    std::cout << "ef.size() = " << ef.size() << '\n';
    std::cout << "ef.back() = " << ef.back() << '\n';
    std::cout << "measured bits/int = " << (8.0 * ef.num_bytes()) / ef.size() << std::endl;

    {
        std::cout << "checking correctness of random access..." << std::endl;
        for (uint64_t i = 0; i != n; ++i) {
            uint64_t got = ef.access(i);  // get the integer at position i
            uint64_t expected = vec[i];
            if (got != expected) {
                std::cout << "got " << got << " at position " << i << "/" << n << " but expected "
                          << expected << std::endl;
                return 1;
            }
        }
        std::cout << "EVERYTHING OK!" << std::endl;
    }

    {
        std::cout << "checking correctness of iterator..." << std::endl;
        auto it = ef.begin();
        uint64_t i = 0;
        while (it.has_next()) {
            uint64_t got = it.next();  // get next integer
            uint64_t expected = vec[i];
            if (got != expected) {
                std::cout << "got " << got << " at position " << i << "/" << n << " but expected "
                          << expected << std::endl;
                return 1;
            }
            ++i;
        }
        assert(i == n);
        std::cout << "EVERYTHING OK!" << std::endl;
    }

    {
        assert(index_zeros == true);
        std::cout << "checking correctness of next_geq..." << std::endl;
        uint64_t i = 0;
        for (auto x : vec) {
            /* since x is in the sequence, next_geq must return (i,x) */
            auto [pos, got] = ef.next_geq(x);
            uint64_t expected = vec[i];

            /* if we have some repeated values, next_neq will return the leftmost occurrence */
            if (i > 0 and vec[i] == vec[i - 1]) {
                ++i;
                continue;
            }

            if (got != expected or pos != i) {
                std::cout << "got " << got << " at position " << pos << "/" << n << " but expected "
                          << expected << " at position " << i << std::endl;
                return 1;
            }
            ++i;
        }
        std::cout << "EVERYTHING OK!" << std::endl;

        std::cout << "checking correctness of next_geq..." << std::endl;
        for (i = 0; i != 10000; ++i) {
            uint64_t x = vec[rand() % n] + (i % 2 == 0 ? 3 : -3);  // get some value
            auto [pos, got] = ef.next_geq(x);
            auto it = std::lower_bound(vec.begin(), vec.end(), x);
            if (it != vec.end()) {
                uint64_t expected = *it;
                uint64_t pos_expected = std::distance(vec.begin(), it);
                if (got != expected or pos != pos_expected) {
                    std::cout << "got " << got << " at position " << pos << "/" << n
                              << " but expected " << expected << " at position " << pos_expected
                              << std::endl;
                    return 1;
                }
            } else {
                assert(x > ef.back());
                uint64_t expected = ef.back();
                if (got != expected or pos != vec.size() - 1) {
                    std::cout << "got " << got << " at position " << pos << "/" << n
                              << " but expected " << expected << " at position " << vec.size() - 1
                              << std::endl;
                    return 1;
                }
            }
        }
        std::cout << "EVERYTHING OK!" << std::endl;
    }

    {
        assert(index_zeros == true);
        std::cout << "checking correctness of prev_leq..." << std::endl;
        uint64_t i = 0;
        for (auto x : vec) {
            /* since x is in the sequence, prev_leq must return i */
            auto pos = ef.prev_leq(x);
            auto got = ef.access(pos);
            uint64_t expected = vec[i];

            /* if we have some repeated values, prev_leq will return the rightmost occurrence */
            if (pos > 0 and vec[pos] == vec[pos - 1]) {
                i = pos + 1;
                continue;
            }

            if (got != expected or pos != i) {
                std::cout << "got " << got << " at position " << pos << "/" << n << " but expected "
                          << expected << " at position " << i << std::endl;
                return 1;
            }
            ++i;
        }
        std::cout << "EVERYTHING OK!" << std::endl;

        std::cout << "checking correctness of prev_leq..." << std::endl;
        const uint64_t front = ef.access(0);
        for (i = 0; i != 10000; ++i) {
            uint64_t x = vec[rand() % n] + (i % 2 == 0 ? 3 : -3);  // get some value
            auto pos = ef.prev_leq(x);
            if (x < front) {
                if (pos != uint64_t(-1)) {
                    std::cout << "expected pos " << uint64_t(-1) << " but got pos = " << pos
                              << std::endl;
                    return 1;
                }
                continue;
            }
            auto got = pos < ef.size() ? ef.access(pos) : ef.back();
            auto it = std::upper_bound(vec.begin(), vec.end(), x) - 1;
            uint64_t expected = *it;
            uint64_t pos_expected = std::distance(vec.begin(), it);
            if (got != expected or pos != pos_expected) {
                std::cout << "got " << got << " at position " << pos << "/" << n << " but expected "
                          << expected << " at position " << pos_expected << std::endl;
                return 1;
            }
        }
        std::cout << "EVERYTHING OK!" << std::endl;
    }

    {
        // test encode_prefix_sum and diff()
        std::cout << "generating other " << n << " integers destributed as Poisson(mean=" << mean
                  << ")..." << std::endl;
        std::generate(vec.begin(), vec.end(), [&]() { return distr(rng); });
        constexpr bool index_zeros = true;
        constexpr bool encode_prefix_sum = true;
        std::cout << "encoding them with elias_fano<index_zeros=" << index_zeros
                  << ", encode_prefix_sum=" << encode_prefix_sum << ">..." << std::endl;
        elias_fano<index_zeros, encode_prefix_sum> ef;
        ef.encode(vec.begin(), vec.size());
        assert(ef.size() == vec.size() + 1);
        std::cout << "ef.size() = " << ef.size() << '\n';
        std::cout << "ef.back() = " << ef.back() << '\n';
        std::cout << "measured bits/int = " << (8.0 * ef.num_bytes()) / ef.size() << std::endl;
        std::cout << "checking correctness of diff..." << std::endl;
        uint64_t i = 0;
        for (auto x : vec) {
            uint64_t got = ef.diff(i);  // get the integer at position i
            uint64_t expected = x;
            if (got != expected) {
                std::cout << "got " << got << " at position " << i << "/" << n << " but expected "
                          << expected << std::endl;
                return 1;
            }
            ++i;
        }
        std::cout << "EVERYTHING OK!" << std::endl;
    }

    return 0;
}