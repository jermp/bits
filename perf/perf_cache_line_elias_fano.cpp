#include "test/common.hpp"
#include "include/cache_line_elias_fano.hpp"
#include "include/elias_fano.hpp"

cache_line_elias_fano encode_with_cache_line_elias_fano(std::vector<uint64_t> const& seq) {
    std::cout << "encoding seq with cache_line_elias_fano..." << std::endl;
    cache_line_elias_fano ef;
    ef.encode(seq.begin(), seq.size());
    REQUIRE(ef.size() == seq.size());
    // std::cout << "ef.size() = " << ef.size() << '\n';
    // std::cout << "ef.back() = " << ef.back() << '\n';
    std::cout << "measured bits/int = " << (8.0 * ef.num_bytes()) / ef.size() << std::endl;
    return ef;
}

elias_fano<false, false> encode_with_elias_fano(std::vector<uint64_t> const& seq) {
    std::cout << "encoding seq with elias_fano<false,false>..." << std::endl;
    elias_fano<false, false> ef;
    ef.encode(seq.begin(), seq.size());
    REQUIRE(ef.size() == seq.size());
    // std::cout << "ef.size() = " << ef.size() << '\n';
    // std::cout << "ef.back() = " << ef.back() << '\n';
    std::cout << "measured bits/int = " << (8.0 * ef.num_bytes()) / ef.size() << std::endl;
    return ef;
}

std::vector<uint64_t> get_queries(const uint64_t num_queries, const uint64_t sequence_length)  //
{
    essentials::uniform_int_rng<uint64_t> distr(0, sequence_length - 1,
                                                essentials::get_random_seed());
    std::vector<uint64_t> queries(num_queries);
    std::generate(queries.begin(), queries.end(), [&]() { return distr.gen(); });
    assert(queries.size() == num_queries);
    return queries;
}

TEST_CASE("access") {
    const int num_runs = 5;
    essentials::timer_type t;

    for (uint64_t log2_sequence_length = 15; log2_sequence_length != 25 + 1;
         ++log2_sequence_length)  //
    {
        uint64_t sequence_length = 1ULL << log2_sequence_length;
        assert(sequence_length > 0);
        std::vector<uint64_t> seq = test::get_sorted_sequence(sequence_length, 300);
        uint64_t num_queries = std::max<uint64_t>(0.1 * sequence_length, 1000);
        std::vector<uint64_t> queries = get_queries(num_queries, sequence_length);
        std::cout << "num_queries = " << num_queries << std::endl;

        {
            t.reset();
            auto ef = encode_with_elias_fano(seq);
            t.start();
            for (int run = 0; run != num_runs; ++run) {
                for (auto i : queries) {
                    uint64_t val = ef.access(i);
                    essentials::do_not_optimize_away(val);
                }
            }
            t.stop();
            std::cout << "  total elapsed time = " << t.elapsed() << std::endl;
            std::cout << "  EF(n=" << sequence_length
                      << ") random access = " << t.elapsed() / (num_runs * num_queries) * 1000.0
                      << " [ns/query]" << std::endl;
        }

        {
            t.reset();
            auto ef = encode_with_cache_line_elias_fano(seq);
            t.start();
            for (int run = 0; run != num_runs; ++run) {
                for (auto i : queries) {
                    uint64_t val = ef.access(i);
                    essentials::do_not_optimize_away(val);
                }
            }
            t.stop();
            std::cout << "  total elapsed time = " << t.elapsed() << std::endl;
            std::cout << "  CLEF(n=" << sequence_length
                      << ") random access = " << t.elapsed() / (num_runs * num_queries) * 1000.0
                      << " [ns/query]" << std::endl;
        }
    }
}