#include "common.hpp"
#include "include/darray.hpp"

constexpr uint64_t num_positions = 10000;

bit_vector encode_with_bit_vector(std::vector<uint64_t> const& seq) {
    std::cout << "encoding seq with a bit_vector..." << std::endl;
    std::cout << "seq.size() = " << seq.size() << std::endl;
    std::cout << "seq.back() = " << seq.back() << std::endl;
    bit_vector::builder bv_builder(seq.back() + 1);
    for (auto x : seq) {
        bv_builder.set(x);          // set bit at position x
        assert(bv_builder.get(x));  // must be set
    }
    bit_vector bv;
    bv_builder.build(bv);
    std::cout << "bv.num_bits() = " << bv.num_bits() << '\n';
    std::cout << "density = " << ((seq.size() + 1) * 100.0) / bv.num_bits() << "%\n";
    std::cout << "measured bits/int = " << (8.0 * bv.num_bytes()) / seq.size() << std::endl;
    return bv;
}

void run_test(const uint64_t max_int) {
    constexpr bool all_distinct = true;
    std::vector<uint64_t> seq = test::get_sorted_sequence(num_positions, max_int, all_distinct);
    auto bv = encode_with_bit_vector(seq);
    darray1 select_index;
    select_index.build(bv);
    for (uint64_t i = 0; i != seq.size(); ++i) {
        uint64_t pos = select_index.select(bv, i);  // position of the i-th set bit
        // std::cout << "select(" << i << ") = " << pos << std::endl;
        REQUIRE_MESSAGE(pos == seq[i], "got " << pos << " but expected " << seq[i]);
    }
}

TEST_CASE("super_sparse") { run_test(32 * 1024); }
TEST_CASE("very_sparse") { run_test(1024); }
TEST_CASE("sparse") { run_test(128); }
TEST_CASE("dense") { run_test(32); }
TEST_CASE("very_dense") { run_test(8); }
TEST_CASE("super_dense") { run_test(2); }
