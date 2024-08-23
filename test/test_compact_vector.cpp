#include "common.hpp"
#include "include/compact_vector.hpp"

constexpr uint64_t sequence_length = 10000;

compact_vector encode_with_compact_vector(std::vector<uint64_t> const& seq) {
    std::cout << "encoding seq with a compact_vector..." << std::endl;
    compact_vector cv;
    cv.build(seq.begin(), seq.size());  // this calls compact_vector::builder::set
    std::cout << "cv.size() = " << cv.size() << '\n';
    std::cout << "cv.width() = " << cv.width() << '\n';
    std::cout << "measured bits/int = " << (8.0 * cv.num_bytes()) / cv.size() << std::endl;
    return cv;
}

TEST_CASE("operator[]") {
    const uint64_t max_int = test::get_random_uint();
    std::cout << "max_int = " << max_int << std::endl;
    std::vector<uint64_t> seq = test::get_sequence(sequence_length, max_int);
    auto cv = encode_with_compact_vector(seq);
    std::cout << "checking correctness of operator[]..." << std::endl;
    for (uint64_t i = 0; i != sequence_length; ++i) {
        uint64_t got = cv[i];  // get the integer at position i
        uint64_t expected = seq[i];
        REQUIRE_MESSAGE(got == expected, "got " << got << " at position " << i << "/"
                                                << sequence_length << " but expected " << expected);
    }
    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("access") {
    const uint64_t max_int = test::get_random_uint();
    std::cout << "max_int = " << max_int << std::endl;
    std::vector<uint64_t> seq = test::get_sequence(sequence_length, max_int);
    auto cv = encode_with_compact_vector(seq);
    std::cout << "checking correctness of access..." << std::endl;
    for (uint64_t i = 0; i != sequence_length; ++i) {
        uint64_t got = cv.access(i);  // get the integer at position i
        uint64_t expected = seq[i];
        REQUIRE_MESSAGE(got == expected, "got " << got << " at position " << i << "/"
                                                << sequence_length << " but expected " << expected);
    }
    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("iterator::operator++") {
    const uint64_t max_int = test::get_random_uint();
    std::cout << "max_int = " << max_int << std::endl;
    std::vector<uint64_t> seq = test::get_sequence(sequence_length, max_int);
    auto cv = encode_with_compact_vector(seq);
    auto cv_it = cv.begin();
    for (uint64_t i = 0; i != seq.size(); ++i, ++cv_it) {
        uint64_t got = *cv_it;
        uint64_t expected = seq[i];
        REQUIRE_MESSAGE(got == expected, "got " << got << " at position " << i << "/" << seq.size()
                                                << " but expected " << expected);
    }
    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("iterator::operator--") {
    const uint64_t max_int = test::get_random_uint();
    std::cout << "max_int = " << max_int << std::endl;
    std::vector<uint64_t> seq = test::get_sequence(sequence_length, max_int);
    auto cv = encode_with_compact_vector(seq);
    auto cv_it = cv.get_iterator_at(sequence_length - 1);
    for (uint64_t i = 0; i != seq.size(); ++i, --cv_it) {
        uint64_t got = *cv_it;
        uint64_t expected = seq[(seq.size() - i) - 1];
        REQUIRE_MESSAGE(got == expected, "got " << got << " at position " << (seq.size() - i) - 1
                                                << "/" << seq.size() << " but expected "
                                                << expected);
    }
    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("iterator::operator+") {
    const uint64_t max_int = test::get_random_uint();
    std::cout << "max_int = " << max_int << std::endl;
    std::vector<uint64_t> seq = test::get_sequence(sequence_length, max_int);
    auto cv = encode_with_compact_vector(seq);
    auto cv_it = cv.begin();
    for (uint64_t i = 0; i != seq.size(); ++i) {
        uint64_t got = *(cv_it + i);
        uint64_t expected = seq[i];
        REQUIRE_MESSAGE(got == expected, "got " << got << " at position " << i << "/" << seq.size()
                                                << " but expected " << expected);
    }
    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("iterator::operator-") {
    const uint64_t max_int = test::get_random_uint();
    std::cout << "max_int = " << max_int << std::endl;
    std::vector<uint64_t> seq = test::get_sequence(sequence_length, max_int);
    auto cv = encode_with_compact_vector(seq);
    auto cv_it = cv.get_iterator_at(sequence_length - 1);
    for (uint64_t i = 0; i != seq.size(); ++i) {
        uint64_t got = *(cv_it - i);
        uint64_t expected = seq[(seq.size() - i) - 1];
        REQUIRE_MESSAGE(got == expected, "got " << got << " at position " << (seq.size() - i) - 1
                                                << "/" << seq.size() << " but expected "
                                                << expected);
    }
    std::cout << "EVERYTHING OK!" << std::endl;
}

TEST_CASE("save_load") {
    const uint64_t max_int = test::get_random_uint();
    std::cout << "max_int = " << max_int << std::endl;
    std::vector<uint64_t> seq = test::get_sequence(sequence_length, max_int);
    auto cv = encode_with_compact_vector(seq);
    const std::string output_filename("cv.bin");
    uint64_t num_saved_bytes = essentials::save(cv, output_filename.c_str());
    std::cout << "num_saved_bytes = " << num_saved_bytes << std::endl;
    compact_vector cv_loaded;
    uint64_t num_loaded_bytes = essentials::load(cv_loaded, output_filename.c_str());
    std::cout << "num_loaded_bytes = " << num_loaded_bytes << std::endl;
    REQUIRE(num_saved_bytes == num_loaded_bytes);
    std::cout << "checking correctness of iterator..." << std::endl;
    auto it = cv_loaded.begin();
    for (uint64_t i = 0; i != seq.size(); ++i, ++it) {
        REQUIRE_MESSAGE(*it == seq[i], "got " << *it << " at position " << i << "/" << seq.size()
                                              << " but expected " << seq[i]);
    }
    std::remove(output_filename.c_str());
    std::cout << "EVERYTHING OK!" << std::endl;
}
