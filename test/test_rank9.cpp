#include "common.hpp"
#include "include/rank9.hpp"

constexpr uint64_t num_positions = 10000;

std::pair<bit_vector, rank9> encode_with_ranked_bit_vector(std::vector<uint64_t> const& seq) {
    std::cout << "encoding seq with a ranked_bit_vector..." << std::endl;
    std::cout << "seq.size() = " << seq.size() << std::endl;
    std::cout << "seq.back() = " << seq.back() << std::endl;
    bit_vector::builder bv_builder(seq.back() + 1);
    for (auto x : seq) {
        bv_builder.set(x);          // set bit at position x
        assert(bv_builder.get(x));  // must be set
    }
    bit_vector bv;
    bv_builder.build(bv);
    rank9 rank_index;
    rank_index.build(bv);
    std::cout << "bv.num_bits() = " << bv.num_bits() << '\n';
    std::cout << "rank_index.num_ones() = " << rank_index.num_ones() << '\n';
    std::cout << "density = " << ((seq.size() + 1) * 100.0) / bv.num_bits() << "%\n";
    std::cout << "bit_vector measured bits/int = " << (8.0 * bv.num_bytes()) / seq.size()
              << std::endl;
    std::cout << "rank_index measured bits/int = " << (8.0 * rank_index.num_bytes()) / seq.size()
              << std::endl;
    return {bv, rank_index};
}

void run_test(const uint64_t max_int) {
    constexpr bool all_distinct = true;
    std::vector<uint64_t> seq = test::get_sorted_sequence(num_positions, max_int, all_distinct);
    auto [B, rank_index] = encode_with_ranked_bit_vector(seq);
    REQUIRE(rank_index.num_ones() == seq.size());
    REQUIRE(rank_index.rank1(B, 0) == 0);
    for (uint64_t i = 0; i != seq.size(); ++i) {
        uint64_t j = seq[i];
        uint64_t num_ones = rank_index.rank1(B, j);   // number of 1s in B[0..j)
        uint64_t num_zeros = rank_index.rank0(B, j);  // number of 0s in B[0..j)
        REQUIRE(num_ones + num_zeros == j);
        // std::cout << "rank1(" << j << ") = " << num_ones << std::endl;
        // std::cout << "rank0(" << j << ") = " << num_zeros << std::endl;
        REQUIRE_MESSAGE(num_ones == i, "got " << num_ones << " but expected " << i);
        REQUIRE_MESSAGE(num_zeros == j - i, "got " << num_ones << " but expected " << j - i);
    }

    const std::string output_filename("r9.bin");
    uint64_t num_saved_bytes = essentials::save(rank_index, output_filename.c_str());
    std::cout << "num_saved_bytes = " << num_saved_bytes << std::endl;
    rank9 rank_index_loaded;
    uint64_t num_loaded_bytes = essentials::load(rank_index_loaded, output_filename.c_str());
    std::cout << "num_loaded_bytes = " << num_loaded_bytes << std::endl;
    std::remove(output_filename.c_str());
    REQUIRE(num_saved_bytes == num_loaded_bytes);

    rank9 other;
    rank_index_loaded.swap(other);
    REQUIRE(other.num_ones() == seq.size());
    REQUIRE(other.rank1(B, 0) == 0);
    for (uint64_t i = 0; i != seq.size(); ++i) {
        uint64_t j = seq[i];
        uint64_t num_ones = other.rank1(B, j);   // number of 1s in B[0..j)
        uint64_t num_zeros = other.rank0(B, j);  // number of 0s in B[0..j)
        REQUIRE(num_ones + num_zeros == j);
        // std::cout << "rank1(" << j << ") = " << num_ones << std::endl;
        // std::cout << "rank0(" << j << ") = " << num_zeros << std::endl;
        REQUIRE_MESSAGE(num_ones == i, "got " << num_ones << " but expected " << i);
        REQUIRE_MESSAGE(num_zeros == j - i, "got " << num_ones << " but expected " << j - i);
    }

    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("super_sparse") { run_test(32 * 1024); }
TEST_CASE("very_sparse") { run_test(1024); }
TEST_CASE("sparse") { run_test(128); }
TEST_CASE("dense") { run_test(32); }
TEST_CASE("very_dense") { run_test(8); }
TEST_CASE("super_dense") { run_test(2); }

void test_save_load_swap(const uint64_t max_int) {
    constexpr bool all_distinct = true;
    std::vector<uint64_t> seq = test::get_sorted_sequence(num_positions, max_int, all_distinct);
    auto [B, rank_index] = encode_with_ranked_bit_vector(seq);

    const std::string output_filename("r9_swap.bin");
    uint64_t num_saved_bytes = 0;

    {
        // Assertions on original
        REQUIRE(rank_index.num_ones() == seq.size());
        REQUIRE(rank_index.rank1(B, 0) == 0);

        num_saved_bytes = essentials::save(rank_index, output_filename.c_str());
        std::cout << "num_saved_bytes = " << num_saved_bytes << std::endl;
    }

    rank9 rank_index_loaded;
    uint64_t num_loaded_bytes = essentials::load(rank_index_loaded, output_filename.c_str());
    std::cout << "num_loaded_bytes = " << num_loaded_bytes << std::endl;
    REQUIRE(num_saved_bytes == num_loaded_bytes);

    std::cout << "checking correctness of rank1 and rank0 after swap..." << std::endl;
    rank9 other;
    rank_index_loaded.swap(other);

    REQUIRE(other.num_ones() == seq.size());
    REQUIRE(other.rank1(B, 0) == 0);
    for (uint64_t i = 0; i != seq.size(); ++i) {
        uint64_t j = seq[i];
        uint64_t num_ones = other.rank1(B, j);   // number of 1s in B[0..j)
        uint64_t num_zeros = other.rank0(B, j);  // number of 0s in B[0..j)
        REQUIRE(num_ones + num_zeros == j);
        REQUIRE_MESSAGE(num_ones == i, "got " << num_ones << " but expected " << i);
        REQUIRE_MESSAGE(num_zeros == j - i, "got " << num_zeros << " but expected " << j - i);
    }

    std::remove(output_filename.c_str());
    std::cout << "EVERYTHING OK!" << std::endl;
}

void test_save_mmap(const uint64_t max_int) {
    constexpr bool all_distinct = true;
    std::vector<uint64_t> seq = test::get_sorted_sequence(num_positions, max_int, all_distinct);
    auto [B, rank_index] = encode_with_ranked_bit_vector(seq);

    const std::string output_filename("r9_mmap.bin");
    uint64_t num_saved_bytes = 0;
    uint64_t num_mapped_bytes = 0;

    {
        num_saved_bytes = essentials::save(rank_index, output_filename.c_str());
        std::cout << "num_saved_bytes = " << num_saved_bytes << std::endl;
    }

    {
        rank9 rank_index_mmapped;
        num_mapped_bytes = essentials::mmap(rank_index_mmapped, output_filename.c_str());
        std::cout << "num_mapped_bytes = " << num_mapped_bytes << std::endl;
        REQUIRE(num_saved_bytes == num_mapped_bytes);

        std::cout << "checking correctness of rank1 and rank0 after mmap..." << std::endl;
        REQUIRE(rank_index_mmapped.num_ones() == seq.size());
        REQUIRE(rank_index_mmapped.rank1(B, 0) == 0);

        for (uint64_t i = 0; i != seq.size(); ++i) {
            uint64_t j = seq[i];
            uint64_t num_ones = rank_index_mmapped.rank1(B, j);
            uint64_t num_zeros = rank_index_mmapped.rank0(B, j);
            REQUIRE(num_ones + num_zeros == j);
            REQUIRE_MESSAGE(num_ones == i, "got " << num_ones << " but expected " << i);
            REQUIRE_MESSAGE(num_zeros == j - i, "got " << num_zeros << " but expected " << j - i);
        }
        std::cout << "EVERYTHING OK!" << std::endl;
    }

    std::remove(output_filename.c_str());
}

TEST_CASE("rank9_save_load_swap_super_sparse") { test_save_load_swap(32 * 1024); }
TEST_CASE("rank9_save_load_swap_very_sparse") { test_save_load_swap(1024); }
TEST_CASE("rank9_save_load_swap_sparse") { test_save_load_swap(128); }
TEST_CASE("rank9_save_load_swap_dense") { test_save_load_swap(32); }
TEST_CASE("rank9_save_load_swap_very_dense") { test_save_load_swap(8); }
TEST_CASE("rank9_save_load_swap_super_dense") { test_save_load_swap(2); }

TEST_CASE("rank9_save_mmap_super_sparse") { test_save_mmap(32 * 1024); }
TEST_CASE("rank9_save_mmap_very_sparse") { test_save_mmap(1024); }
TEST_CASE("rank9_save_mmap_sparse") { test_save_mmap(128); }
TEST_CASE("rank9_save_mmap_dense") { test_save_mmap(32); }
TEST_CASE("rank9_save_mmap_very_dense") { test_save_mmap(8); }
TEST_CASE("rank9_save_mmap_super_dense") { test_save_mmap(2); }
