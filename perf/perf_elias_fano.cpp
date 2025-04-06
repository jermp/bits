#include "test/common.hpp"
#include "perf/common.hpp"
#include "include/elias_fano.hpp"

#include <sstream>

elias_fano<false, false> encode_with_elias_fano(std::vector<uint64_t> const& seq) {
    std::cout << "encoding seq with elias_fano<false,false>..." << std::endl;
    elias_fano<false, false> ef;
    ef.encode(seq.begin(), seq.size());
    REQUIRE(ef.size() == seq.size());
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
        auto ef = encode_with_elias_fano(seq);
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