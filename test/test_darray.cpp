#include "common.hpp"
#include "include/darray.hpp"

constexpr uint64_t num_positions = 10000;

bit_vector encode_with_bit_vector(std::vector<uint64_t> const& seq, bool index_ones) {
    std::cout << "encoding seq with a bit_vector..." << std::endl;
    std::cout << "seq.size() = " << seq.size() << std::endl;
    std::cout << "seq.back() = " << seq.back() << std::endl;
    // bit_vector::builder bv_builder(50); // example from notes
    bit_vector::builder bv_builder(seq.back() + 1, !index_ones);
    for (auto x : seq) {
        bv_builder.set(x, index_ones);            // set bit at position x
        assert(bv_builder.get(x) == index_ones);  // must be set
    }
    bit_vector bv;
    bv_builder.build(bv);
    std::cout << "bv.num_bits() = " << bv.num_bits() << '\n';
    std::cout << "density = " << ((seq.size() + 1) * 100.0) / bv.num_bits() << "%\n";
    std::cout << "measured bits/int = " << (8.0 * bv.num_bytes()) / seq.size() << std::endl;
    return bv;
}

template <typename DArray, bool index_ones>
void run_test(const uint64_t max_int) {
    constexpr bool all_distinct = true;
    // Example from notes.
    // std::vector<uint64_t> seq = {0,  1,  3,  5,  12, 15, 18, 19, 20,
    //                              24, 29, 34, 36, 39, 40, 41, 42, 47};
    std::vector<uint64_t> seq = test::get_sorted_sequence(num_positions, max_int, all_distinct);
    auto bv = encode_with_bit_vector(seq, index_ones);
    DArray select_index;
    select_index.build(bv);
    REQUIRE(select_index.num_positions() == seq.size());
    for (uint64_t i = 0; i != seq.size(); ++i) {
        uint64_t pos = select_index.select(bv, i);  // position of the i-th bit set
        // std::cout << "select(" << i << ") = " << pos << std::endl;
        REQUIRE_MESSAGE(pos == seq[i], "got " << pos << " but expected " << seq[i]);
    }

    const std::string output_filename("darray.bin");
    uint64_t num_saved_bytes = essentials::save(select_index, output_filename.c_str());
    std::cout << "num_saved_bytes = " << num_saved_bytes << std::endl;
    DArray select_index_loaded;
    uint64_t num_loaded_bytes = essentials::load(select_index_loaded, output_filename.c_str());
    std::cout << "num_loaded_bytes = " << num_loaded_bytes << std::endl;
    std::remove(output_filename.c_str());
    REQUIRE(num_saved_bytes == num_loaded_bytes);

    DArray other;
    select_index_loaded.swap(other);
    REQUIRE(select_index.num_positions() == seq.size());
    for (uint64_t i = 0; i != seq.size(); ++i) {
        uint64_t pos = other.select(bv, i);  // position of the i-th bit set
        // std::cout << "select(" << i << ") = " << pos << std::endl;
        REQUIRE_MESSAGE(pos == seq[i], "got " << pos << " but expected " << seq[i]);
    }

    std::cout << "EVERYTHING OK!" << std::endl;
}

// TEST_CASE("example-from-notes") { run_test<darray<util::identity_getter, 4, 2>, true>(0); }

TEST_CASE("super_sparse1") { run_test<darray1, true>(32 * 1024); }
TEST_CASE("super_sparse0") { run_test<darray0, false>(32 * 1024); }

TEST_CASE("very_sparse1") { run_test<darray1, true>(1024); }
TEST_CASE("very_sparse0") { run_test<darray0, false>(1024); }

TEST_CASE("sparse1") { run_test<darray1, true>(128); }
TEST_CASE("sparse0") { run_test<darray0, false>(128); }

TEST_CASE("dense1") { run_test<darray1, true>(32); }
TEST_CASE("dense0") { run_test<darray0, false>(32); }

TEST_CASE("very_dense1") { run_test<darray1, true>(8); }
TEST_CASE("very_dense0") { run_test<darray0, false>(8); }

TEST_CASE("super_dense1") { run_test<darray1, true>(2); }
TEST_CASE("super_dense0") { run_test<darray0, false>(2); }

template <typename DArray, bool index_ones>
void test_save_load_swap(const uint64_t max_int) {
    constexpr bool all_distinct = true;
    std::vector<uint64_t> seq = test::get_sorted_sequence(num_positions, max_int, all_distinct);
    auto bv = encode_with_bit_vector(seq, index_ones);

    const std::string output_filename("darray_swap.bin");
    uint64_t num_saved_bytes = 0;

    {
        DArray select_index;
        select_index.build(bv);
        REQUIRE(select_index.num_positions() == seq.size());

        // Save to disk
        num_saved_bytes = essentials::save(select_index, output_filename.c_str());
        std::cout << "num_saved_bytes = " << num_saved_bytes << std::endl;
    }

    DArray select_index_loaded;
    uint64_t num_loaded_bytes = essentials::load(select_index_loaded, output_filename.c_str());
    std::cout << "num_loaded_bytes = " << num_loaded_bytes << std::endl;
    REQUIRE(num_saved_bytes == num_loaded_bytes);

    std::cout << "checking correctness of swap and select..." << std::endl;
    DArray other;
    select_index_loaded.swap(other);

    REQUIRE(other.num_positions() == seq.size());
    for (uint64_t i = 0; i != seq.size(); ++i) {
        uint64_t pos = other.select(bv, i);
        REQUIRE_MESSAGE(pos == seq[i], "got " << pos << " but expected " << seq[i]);
    }

    std::remove(output_filename.c_str());
    std::cout << "EVERYTHING OK!" << std::endl;
}

template <typename DArray, bool index_ones>
void test_save_mmap(const uint64_t max_int) {
    constexpr bool all_distinct = true;
    std::vector<uint64_t> seq = test::get_sorted_sequence(num_positions, max_int, all_distinct);

    // The bit_vector must remain alive for the mmapped select_index to use it
    auto bv = encode_with_bit_vector(seq, index_ones);

    const std::string output_filename("darray_mmap.bin");
    uint64_t num_saved_bytes = 0;
    uint64_t num_mapped_bytes = 0;

    {
        DArray select_index;
        select_index.build(bv);
        REQUIRE(select_index.num_positions() == seq.size());

        num_saved_bytes = essentials::save(select_index, output_filename.c_str());
        std::cout << "num_saved_bytes = " << num_saved_bytes << std::endl;
    }

    {
        DArray select_index_mmapped;
        num_mapped_bytes = essentials::mmap(select_index_mmapped, output_filename.c_str());
        std::cout << "num_mapped_bytes = " << num_mapped_bytes << std::endl;
        REQUIRE(num_saved_bytes == num_mapped_bytes);

        std::cout << "checking correctness of mmapped select..." << std::endl;
        REQUIRE(select_index_mmapped.num_positions() == seq.size());

        for (uint64_t i = 0; i != seq.size(); ++i) {
            uint64_t pos = select_index_mmapped.select(bv, i);
            REQUIRE_MESSAGE(pos == seq[i], "got " << pos << " but expected " << seq[i]);
        }
        std::cout << "EVERYTHING OK!" << std::endl;
    }

    std::remove(output_filename.c_str());
}

TEST_CASE("darray_save_load_swap_super_sparse1") { test_save_load_swap<darray1, true>(32 * 1024); }
TEST_CASE("darray_save_load_swap_super_sparse0") { test_save_load_swap<darray0, false>(32 * 1024); }

TEST_CASE("darray_save_load_swap_dense1") { test_save_load_swap<darray1, true>(32); }
TEST_CASE("darray_save_load_swap_dense0") { test_save_load_swap<darray0, false>(32); }

TEST_CASE("darray_save_mmap_super_sparse1") { test_save_mmap<darray1, true>(32 * 1024); }
TEST_CASE("darray_save_mmap_super_sparse0") { test_save_mmap<darray0, false>(32 * 1024); }

TEST_CASE("darray_save_mmap_dense1") { test_save_mmap<darray1, true>(32); }
TEST_CASE("darray_save_mmap_dense0") { test_save_mmap<darray0, false>(32); }
