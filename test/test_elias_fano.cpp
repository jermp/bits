#include "common.hpp"
#include "include/elias_fano.hpp"

constexpr uint64_t sequence_length = 10000;

template <bool index_zeros, bool encode_prefix_sum>
elias_fano<index_zeros, encode_prefix_sum> encode_with_elias_fano(
    std::vector<uint64_t> const& seq) {
    std::cout << "encoding them with elias_fano<index_zeros=" << index_zeros
              << ", encode_prefix_sum=" << encode_prefix_sum << ">..." << std::endl;
    elias_fano<index_zeros, encode_prefix_sum> ef;
    ef.encode(seq.begin(), seq.size());
    std::cout << "ef.size() = " << ef.size() << '\n';
    std::cout << "ef.back() = " << ef.back() << '\n';
    std::cout << "measured bits/int = " << (8.0 * ef.num_bytes()) / ef.size() << std::endl;
    return ef;
}

TEST_CASE("access") {
    std::vector<uint64_t> seq = bits::test::get_sorted_sequence(sequence_length);
    constexpr bool index_zeros = true;
    constexpr bool encode_prefix_sum = false;
    auto ef = encode_with_elias_fano<index_zeros, encode_prefix_sum>(seq);
    std::cout << "checking correctness of random access..." << std::endl;
    for (uint64_t i = 0; i != sequence_length; ++i) {
        uint64_t got = ef.access(i);  // get the integer at position i
        uint64_t expected = seq[i];
        REQUIRE_MESSAGE(got == expected, "got " << got << " at position " << i << "/"
                                                << sequence_length << " but expected " << expected);
    }
    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("iterator") {
    std::vector<uint64_t> seq = bits::test::get_sorted_sequence(sequence_length);
    constexpr bool index_zeros = true;
    constexpr bool encode_prefix_sum = false;
    auto ef = encode_with_elias_fano<index_zeros, encode_prefix_sum>(seq);
    std::cout << "checking correctness of iterator..." << std::endl;
    auto it = ef.begin();
    uint64_t i = 0;
    while (it.has_next()) {
        uint64_t got = it.next();  // get next integer
        uint64_t expected = seq[i];
        REQUIRE_MESSAGE(got == expected, "got " << got << " at position " << i << "/"
                                                << sequence_length << " but expected " << expected);
        ++i;
    }
    assert(i == sequence_length);
    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("next_geq") {
    std::vector<uint64_t> seq = bits::test::get_sorted_sequence(sequence_length);
    constexpr bool index_zeros = true;
    constexpr bool encode_prefix_sum = false;
    auto ef = encode_with_elias_fano<index_zeros, encode_prefix_sum>(seq);
    std::cout << "checking correctness of next_geq..." << std::endl;
    uint64_t i = 0;
    for (auto x : seq) {
        /* since x is in the sequence, next_geq must return (i,x) */

        // auto [pos, got] = ef.next_geq(x); // can't use structural binding with doctest...
        std::pair<uint64_t, uint64_t> p = ef.next_geq(x);
        uint64_t pos = p.first;
        uint64_t got = p.second;
        uint64_t expected = seq[i];

        // std::cout << "x=" << x << "; pos=" << pos << "; got=" << got << std::endl;

        /* if we have some repeated values, next_neq will return the leftmost occurrence */
        if (i > 0 and seq[i] == seq[i - 1]) {
            ++i;
            continue;
        }

        bool good = (got == expected) && (pos == i);
        REQUIRE_MESSAGE(good, "got " << got << " at position " << pos << "/" << sequence_length
                                     << " but expected " << expected << " at position " << i);
        ++i;
    }
    std::cout << "EVERYTHING OK!" << std::endl;

    std::cout << "checking correctness of next_geq..." << std::endl;
    for (i = 0; i != 10000; ++i) {
        uint64_t x = seq[rand() % sequence_length] + (i % 2 == 0 ? 3 : -3);  // get some value

        // auto [pos, got] = ef.next_geq(x); // can't use structural binding with doctest...
        std::pair<uint64_t, uint64_t> p = ef.next_geq(x);
        uint64_t pos = p.first;
        uint64_t got = p.second;

        // std::cout << "x=" << x << "; pos=" << pos << "; got=" << got << std::endl;
        auto it = std::lower_bound(seq.begin(), seq.end(), x);
        if (it != seq.end()) {
            uint64_t expected = *it;
            uint64_t pos_expected = std::distance(seq.begin(), it);
            bool good = (got == expected) && (pos == pos_expected);
            REQUIRE_MESSAGE(good, "got " << got << " at position " << pos << "/" << sequence_length
                                         << " but expected " << expected << " at position "
                                         << pos_expected);
        } else {
            assert(x > ef.back());
            uint64_t expected = ef.back();
            bool good = (got == expected) && (pos == seq.size() - 1);
            REQUIRE_MESSAGE(good, "got " << got << " at position " << pos << "/" << sequence_length
                                         << " but expected " << expected << " at position "
                                         << seq.size() - 1);
        }
    }
    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("small_next_geq") {
    std::vector<uint64_t> seq{1, 1, 1, 3, 3, 5, 6, 6, 6, 9, 13, 23, 23, 23};
    constexpr bool index_zeros = true;
    constexpr bool encode_prefix_sum = false;
    auto ef = encode_with_elias_fano<index_zeros, encode_prefix_sum>(seq);
    std::pair<uint64_t, uint64_t> p;
    uint64_t pos, got;

    p = ef.next_geq(1);
    pos = p.first;
    got = p.second;
    REQUIRE(pos == 0);  // leftmost
    REQUIRE(got == 1);

    p = ef.next_geq(3);
    pos = p.first;
    got = p.second;
    REQUIRE(pos == 3);  // leftmost
    REQUIRE(got == 3);

    p = ef.next_geq(6);
    pos = p.first;
    got = p.second;
    REQUIRE(pos == 6);  // leftmost
    REQUIRE(got == 6);

    p = ef.next_geq(23);
    pos = p.first;
    got = p.second;
    REQUIRE(pos == 11);  // leftmost
    REQUIRE(got == 23);

    p = ef.next_geq(100);  // saturate
    pos = p.first;
    got = p.second;
    REQUIRE(pos == seq.size() - 1);
    REQUIRE(got == seq.back());
}

TEST_CASE("small_prev_leq") {
    std::vector<uint64_t> seq{1, 1, 1, 3, 3, 5, 6, 6, 6, 9, 13, 23, 23, 23};
    constexpr bool index_zeros = true;
    constexpr bool encode_prefix_sum = false;
    auto ef = encode_with_elias_fano<index_zeros, encode_prefix_sum>(seq);
    uint64_t pos;

    pos = ef.prev_leq(0);
    REQUIRE(pos == uint64_t(-1));  // undefined since 0 < 1

    pos = ef.prev_leq(1);
    REQUIRE(pos == 2);  // rightmost

    pos = ef.prev_leq(3);
    REQUIRE(pos == 4);  // rightmost

    pos = ef.prev_leq(6);
    REQUIRE(pos == 8);  // rightmost

    pos = ef.prev_leq(23);
    REQUIRE(pos == seq.size() - 1);  // rightmost

    pos = ef.prev_leq(100);  // saturate
    REQUIRE(pos == seq.size() - 1);
}

TEST_CASE("prev_leq") {
    std::vector<uint64_t> seq = bits::test::get_sorted_sequence(sequence_length);
    constexpr bool index_zeros = true;
    constexpr bool encode_prefix_sum = false;
    auto ef = encode_with_elias_fano<index_zeros, encode_prefix_sum>(seq);
    std::cout << "checking correctness of prev_leq..." << std::endl;
    uint64_t i = 0;
    for (auto x : seq) {
        /* since x is in the sequence, prev_leq must return i */
        auto pos = ef.prev_leq(x);
        auto got = ef.access(pos);
        uint64_t expected = seq[i];

        // std::cout << "x=" << x << "; pos=" << pos << "; got=" << got << std::endl;

        /* if we have some repeated values, prev_leq will return the rightmost occurrence */
        if (pos > 0 and seq[pos] == seq[pos - 1]) {
            i = pos + 1;
            continue;
        }
        bool good = (got == expected) && (pos == i);
        REQUIRE_MESSAGE(good, "got " << got << " at position " << pos << "/" << sequence_length
                                     << " but expected " << expected << " at position " << i);
        ++i;
    }
    std::cout << "EVERYTHING OK!" << std::endl;

    std::cout << "checking correctness of prev_leq..." << std::endl;
    const uint64_t front = ef.access(0);
    for (i = 0; i != 10000; ++i) {
        uint64_t x = seq[rand() % sequence_length] + (i % 2 == 0 ? 3 : -3);  // get some value
        auto pos = ef.prev_leq(x);
        if (x < front) {
            REQUIRE_MESSAGE(pos == uint64_t(-1),
                            "expected pos " << uint64_t(-1) << " but got pos = " << pos);
            continue;
        }
        auto got = pos < ef.size() ? ef.access(pos) : ef.back();
        // std::cout << "x=" << x << "; pos=" << pos << "; got=" << got << std::endl;
        auto it = std::upper_bound(seq.begin(), seq.end(), x) - 1;
        uint64_t expected = *it;
        uint64_t pos_expected = std::distance(seq.begin(), it);
        bool good = (got == expected) && (pos == pos_expected);
        REQUIRE_MESSAGE(good, "got " << got << " at position " << pos << "/" << sequence_length
                                     << " but expected " << expected << " at position "
                                     << pos_expected);
    }
    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("diff") {
    std::vector<uint64_t> seq = bits::test::get_sequence(sequence_length);
    constexpr bool index_zeros = true;
    constexpr bool encode_prefix_sum = true;
    auto ef = encode_with_elias_fano<index_zeros, encode_prefix_sum>(seq);
    assert(ef.size() == seq.size() + 1);
    std::cout << "ef.size() = " << ef.size() << '\n';
    std::cout << "ef.back() = " << ef.back() << '\n';
    std::cout << "measured bits/int = " << (8.0 * ef.num_bytes()) / ef.size() << std::endl;
    std::cout << "checking correctness of diff..." << std::endl;
    uint64_t i = 0;
    for (auto x : seq) {
        uint64_t got = ef.diff(i);  // get the integer at position i
        uint64_t expected = x;
        REQUIRE_MESSAGE(got == expected, "got " << got << " at position " << i << "/"
                                                << sequence_length << " but expected " << expected);
        ++i;
    }
    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("save_load") {
    std::vector<uint64_t> seq = bits::test::get_sorted_sequence(sequence_length);
    constexpr bool index_zeros = true;
    constexpr bool encode_prefix_sum = false;
    using ef_type = elias_fano<index_zeros, encode_prefix_sum>;
    ef_type ef = encode_with_elias_fano<index_zeros, encode_prefix_sum>(seq);
    const std::string output_filename("ef.bin");
    uint64_t num_saved_bytes = essentials::save(ef, output_filename.c_str());
    std::cout << "num_saved_bytes = " << num_saved_bytes << std::endl;
    ef_type ef_loaded;
    uint64_t num_loaded_bytes = essentials::load(ef_loaded, output_filename.c_str());
    std::cout << "num_loaded_bytes = " << num_loaded_bytes << std::endl;
    REQUIRE(num_saved_bytes == num_loaded_bytes);
    std::cout << "checking correctness of iterator..." << std::endl;
    auto it = ef_loaded.begin();
    uint64_t i = 0;
    while (it.has_next()) {
        uint64_t got = it.next();  // get next integer
        uint64_t expected = seq[i];
        REQUIRE_MESSAGE(got == expected, "got " << got << " at position " << i << "/"
                                                << sequence_length << " but expected " << expected);
        ++i;
    }
    assert(i == sequence_length);
    std::cout << "EVERYTHING OK!" << std::endl;
    std::remove(output_filename.c_str());
}
