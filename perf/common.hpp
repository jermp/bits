#pragma once

namespace bits::perf {

const int num_runs = 5;
const uint64_t seed = 0;
const uint64_t min_log2_sequence_length = 15;
const uint64_t max_log2_sequence_length = 25;

std::vector<uint64_t> get_queries(const uint64_t num_queries, const uint64_t sequence_length,
                                  const uint64_t seed = essentials::get_random_seed())  //
{
    essentials::uniform_int_rng<uint64_t> distr(0, sequence_length - 1, seed);
    std::vector<uint64_t> queries(num_queries);
    std::generate(queries.begin(), queries.end(), [&]() { return distr.gen(); });
    assert(queries.size() == num_queries);
    return queries;
}

}  // namespace bits::perf