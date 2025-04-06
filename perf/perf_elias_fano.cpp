#include "test/common.hpp"
#include "perf/common.hpp"
#include "include/elias_fano.hpp"

#include <sstream>

template <bool index_zeros, bool encode_prefix_sum>
elias_fano<index_zeros, encode_prefix_sum>                //
encode_with_elias_fano(std::vector<uint64_t> const& seq)  //
{
    std::cout << "encoding seq with elias_fano<false,false>..." << std::endl;
    elias_fano<index_zeros, encode_prefix_sum> ef;
    ef.encode(seq.begin(), seq.size());
    REQUIRE(ef.size() == seq.size() + encode_prefix_sum);
    std::cout << "measured bits/int = " << (8.0 * ef.num_bytes()) / ef.size() << std::endl;
    return ef;
}

TEST_CASE("access") {
    essentials::timer_type t;

    std::stringstream ss_avg_ns_per_query;
    std::stringstream ss_sequence_lengths;
    ss_avg_ns_per_query << "\"avg_ns_per_query\":[";
    ss_sequence_lengths << "\"sequence_lengths\":[";

    for (uint64_t log2_sequence_length = perf::min_log2_sequence_length;
         log2_sequence_length <= perf::max_log2_sequence_length; ++log2_sequence_length)  //
    {
        uint64_t sequence_length = 1ULL << log2_sequence_length;
        assert(sequence_length > 0);
        constexpr bool all_distinct = false;
        std::vector<uint64_t> seq =
            test::get_sorted_sequence(sequence_length, 300, all_distinct, perf::seed);
        // test::print(seq);
        uint64_t num_queries = std::max<uint64_t>(0.1 * sequence_length, 1000);
        std::vector<uint64_t> queries = perf::get_queries(num_queries, sequence_length, perf::seed);
        // test::print(queries);
        t.reset();
        auto ef = encode_with_elias_fano<false, false>(seq);
        t.start();
        for (int run = 0; run != perf::num_runs; ++run) {
            for (auto i : queries) {
                uint64_t val = ef.access(i);
                essentials::do_not_optimize_away(val);
            }
        }
        t.stop();
        double avg_ns_per_query = t.elapsed() / (perf::num_runs * num_queries) * 1000.0;
        std::cout << "  total elapsed time = " << t.elapsed() << std::endl;
        std::cout << "  EF(n=" << sequence_length << ") random access = " << avg_ns_per_query
                  << " [ns/query]" << std::endl;

        ss_sequence_lengths << sequence_length;
        ss_avg_ns_per_query << avg_ns_per_query;
        if (log2_sequence_length != perf::max_log2_sequence_length) {
            ss_avg_ns_per_query << ',';
            ss_sequence_lengths << ',';
        }
    }

    ss_avg_ns_per_query << ']';
    ss_sequence_lengths << ']';
    std::string json("{\"query\":\"elias_fano::access\", \"seed\":" + std::to_string(perf::seed) +
                     ", ");
    json += ss_sequence_lengths.str();
    json.push_back(',');
    json.push_back(' ');
    json += ss_avg_ns_per_query.str();
    json.push_back('}');
    std::cerr << json << std::endl;
}

TEST_CASE("diff") {
    essentials::timer_type t;

    std::stringstream ss_avg_ns_per_query;
    std::stringstream ss_sequence_lengths;
    ss_avg_ns_per_query << "\"avg_ns_per_query\":[";
    ss_sequence_lengths << "\"sequence_lengths\":[";

    for (uint64_t log2_sequence_length = perf::min_log2_sequence_length;
         log2_sequence_length <= perf::max_log2_sequence_length; ++log2_sequence_length)  //
    {
        uint64_t sequence_length = 1ULL << log2_sequence_length;
        assert(sequence_length > 0);
        std::vector<uint64_t> seq = test::get_sequence(sequence_length, 300, perf::seed);
        // test::print(seq);
        uint64_t num_queries = std::max<uint64_t>(0.1 * sequence_length, 1000);
        std::vector<uint64_t> queries = perf::get_queries(num_queries, sequence_length, perf::seed);
        // test::print(queries);
        t.reset();
        auto ef = encode_with_elias_fano<false, true>(seq);
        t.start();
        for (int run = 0; run != perf::num_runs; ++run) {
            for (auto i : queries) {
                uint64_t val = ef.diff(i);
                essentials::do_not_optimize_away(val);
            }
        }
        t.stop();
        double avg_ns_per_query = t.elapsed() / (perf::num_runs * num_queries) * 1000.0;
        std::cout << "  total elapsed time = " << t.elapsed() << std::endl;
        std::cout << "  EF(n=" << sequence_length << ") diff = " << avg_ns_per_query
                  << " [ns/query]" << std::endl;

        ss_sequence_lengths << sequence_length;
        ss_avg_ns_per_query << avg_ns_per_query;
        if (log2_sequence_length != perf::max_log2_sequence_length) {
            ss_avg_ns_per_query << ',';
            ss_sequence_lengths << ',';
        }
    }

    ss_avg_ns_per_query << ']';
    ss_sequence_lengths << ']';
    std::string json("{\"query\":\"elias_fano::diff\", \"seed\":" + std::to_string(perf::seed) +
                     ", ");
    json += ss_sequence_lengths.str();
    json.push_back(',');
    json.push_back(' ');
    json += ss_avg_ns_per_query.str();
    json.push_back('}');
    std::cerr << json << std::endl;
}

template <typename Func>
void perf_func(std::string const& what, double const avg_gap, Func f) {
    essentials::timer_type t;

    std::stringstream ss_avg_ns_per_query;
    std::stringstream ss_sequence_lengths;
    ss_avg_ns_per_query << "\"avg_ns_per_query\":[";
    ss_sequence_lengths << "\"sequence_lengths\":[";

    for (uint64_t log2_sequence_length = perf::min_log2_sequence_length;
         log2_sequence_length <= perf::max_log2_sequence_length; ++log2_sequence_length)  //
    {
        uint64_t sequence_length = 1ULL << log2_sequence_length;
        assert(sequence_length > 0);
        const uint64_t universe = avg_gap * sequence_length;
        std::vector<uint64_t> seq =
            test::get_uniform_sorted_sequence(sequence_length, universe, perf::seed);
        // test::print(seq);
        uint64_t num_queries = std::max<uint64_t>(0.1 * sequence_length, 1000);
        std::vector<uint64_t> queries = perf::get_queries(num_queries, universe + 1, perf::seed);
        // test::print(queries);
        t.reset();
        auto ef = encode_with_elias_fano<true, false>(seq);
        t.start();
        for (int run = 0; run != perf::num_runs; ++run) {
            for (auto x : queries) {
                auto ret = f(ef, x);
                essentials::do_not_optimize_away(&ret);
            }
        }
        t.stop();
        double avg_ns_per_query = t.elapsed() / (perf::num_runs * num_queries) * 1000.0;
        std::cout << "  total elapsed time = " << t.elapsed() << std::endl;
        std::cout << "  EF(n=" << sequence_length << ") next_geq = " << avg_ns_per_query
                  << " [ns/query]" << std::endl;

        ss_sequence_lengths << sequence_length;
        ss_avg_ns_per_query << avg_ns_per_query;
        if (log2_sequence_length != perf::max_log2_sequence_length) {
            ss_avg_ns_per_query << ',';
            ss_sequence_lengths << ',';
        }
    }

    ss_avg_ns_per_query << ']';
    ss_sequence_lengths << ']';
    std::string json("{\"query\":\"elias_fano::" + what +
                     "\", \"seed\":" + std::to_string(perf::seed) + ", ");
    json += ss_sequence_lengths.str();
    json.push_back(',');
    json.push_back(' ');
    json += ss_avg_ns_per_query.str();
    json.push_back('}');
    std::cerr << json << std::endl;
}

TEST_CASE("next_geq_dense")  //
{
    perf_func("next_geq_dense", 3.0,
              [](elias_fano<true, false> const& ef, const uint64_t x) { return ef.next_geq(x); });
}

TEST_CASE("next_geq_sparse")  //
{
    perf_func("next_geq_sparse", 3000.0,
              [](elias_fano<true, false> const& ef, const uint64_t x) { return ef.next_geq(x); });
}

TEST_CASE("prev_leq_dense")  //
{
    perf_func("prev_leq_dense", 3.0,
              [](elias_fano<true, false> const& ef, const uint64_t x) { return ef.prev_leq(x); });
}

TEST_CASE("prev_leq_sparse")  //
{
    perf_func("prev_leq_sparse", 3000.0,
              [](elias_fano<true, false> const& ef, const uint64_t x) { return ef.prev_leq(x); });
}

TEST_CASE("locate_dense")  //
{
    perf_func("locate_dense", 3.0,
              [](elias_fano<true, false> const& ef, const uint64_t x) { return ef.locate(x); });
}

TEST_CASE("locate_sparse")  //
{
    perf_func("locate_sparse", 3000.0,
              [](elias_fano<true, false> const& ef, const uint64_t x) { return ef.locate(x); });
}
