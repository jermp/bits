#include "common.hpp"
#include "include/ranked_bit_vector.hpp"

constexpr uint64_t num_positions = 10000;

ranked_bit_vector encode_with_ranked_bit_vector(std::vector<uint64_t> const& seq) {
    std::cout << "encoding seq with a ranked_bit_vector..." << std::endl;
    std::cout << "seq.size() = " << seq.size() << std::endl;
    std::cout << "seq.back() = " << seq.back() << std::endl;
    ranked_bit_vector::builder rbv_builder(seq.back() + 1);
    for (auto x : seq) {
        rbv_builder.set(x);          // set bit at position x
        assert(rbv_builder.get(x));  // must be set
    }
    ranked_bit_vector rbv;
    rbv_builder.build(rbv);
    rbv.build_index();
    std::cout << "rbv.num_bits() = " << rbv.num_bits() << '\n';
    std::cout << "rbv.num_ones() = " << rbv.num_ones() << '\n';
    std::cout << "rbv.num_zeros() = " << rbv.num_zeros() << '\n';
    std::cout << "density = " << ((seq.size() + 1) * 100.0) / rbv.num_bits() << "%\n";
    std::cout << "measured bits/int = " << (8.0 * rbv.num_bytes()) / seq.size() << std::endl;
    return rbv;
}

void run_test(const uint64_t max_int) {
    constexpr bool all_distinct = true;
    std::vector<uint64_t> seq = test::get_sorted_sequence(num_positions, max_int, all_distinct);
    auto rbv = encode_with_ranked_bit_vector(seq);
    REQUIRE(rbv.num_ones() == seq.size());
    REQUIRE(rbv.rank1(0) == 0);
    for (uint64_t i = 0; i != seq.size(); ++i) {
        uint64_t j = seq[i];
        uint64_t num_ones = rbv.rank1(j);   // number of 1s in B[0..j)
        uint64_t num_zeros = rbv.rank0(j);  // number of 0s in B[0..j)
        REQUIRE(num_ones + num_zeros == j);
        // std::cout << "rank1(" << j << ") = " << num_ones << std::endl;
        // std::cout << "rank0(" << j << ") = " << num_zeros << std::endl;
        REQUIRE_MESSAGE(num_ones == i, "got " << num_ones << " but expected " << i);
        REQUIRE_MESSAGE(num_zeros == j - i, "got " << num_ones << " but expected " << j - i);
    }
}

TEST_CASE("super_sparse") { run_test(32 * 1024); }
TEST_CASE("very_sparse") { run_test(1024); }
TEST_CASE("sparse") { run_test(128); }
TEST_CASE("dense") { run_test(32); }
TEST_CASE("very_dense") { run_test(8); }
TEST_CASE("super_dense") { run_test(2); }
