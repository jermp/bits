#include "common.hpp"
#include "include/integer_codes.hpp"

constexpr uint64_t sequence_length = 10000;

TEST_CASE("32bits") {
    std::vector<uint64_t> seq = bits::test::get_sequence(sequence_length);
    bit_vector::builder builder;
    constexpr uint64_t max_uint32 = uint64_t(1) << 32;
    for (auto x : seq) {
        if (x < max_uint32) util::write_32bits(builder, x);
    }
    std::cout << "measured bits/int = " << (1.0 * builder.num_bits()) / seq.size() << std::endl;

    // directly read the bits of the bit_vector::builder without building a bit_vector
    bit_vector_iterator bv_it(builder.data().data(), builder.data().size());

    // or first build a bit_vector...
    // bit_vector bv;
    // builder.build(bv);
    // bit_vector_iterator bv_it(bv.data().data(), bv.data().size());

    for (auto x : seq) {
        if (x < max_uint32) {
            uint64_t got = util::read_32bits(bv_it);
            REQUIRE_MESSAGE(got == x, "got " << got << " but expected " << x);
        }
    }
}

TEST_CASE("unary") {
    std::vector<uint64_t> seq =
        bits::test::get_sequence(sequence_length, 64 - 1  // each value must be < 64
        );
    bit_vector::builder builder;
    for (auto x : seq) util::write_unary(builder, x);
    std::cout << "measured bits/int = " << (1.0 * builder.num_bits()) / seq.size() << std::endl;
    bit_vector_iterator bv_it(builder.data().data(), builder.data().size());
    for (auto x : seq) {
        uint64_t got = util::read_unary(bv_it);
        REQUIRE_MESSAGE(got == x, "got " << got << " but expected " << x);
    }
}

TEST_CASE("binary") {
    const uint64_t u = rand() % 1000 + 1;  // at least 1
    std::cout << "u = " << u << std::endl;
    std::vector<uint64_t> seq = bits::test::get_sequence(sequence_length, u);
    bit_vector::builder builder;
    for (auto x : seq) {
        assert(x <= u);
        util::write_binary(builder, x, u);
    }
    std::cout << "measured bits/int = " << (1.0 * builder.num_bits()) / seq.size() << std::endl;
    bit_vector_iterator bv_it(builder.data().data(), builder.data().size());
    for (auto x : seq) {
        uint64_t got = util::read_binary(bv_it, u);
        REQUIRE_MESSAGE(got == x, "got " << got << " but expected " << x);
    }
}

TEST_CASE("gamma") {
    const uint64_t u = rand() % 1000 + 1;  // at least 1
    std::cout << "u = " << u << std::endl;
    std::vector<uint64_t> seq = bits::test::get_sequence(sequence_length, u);
    bit_vector::builder builder;
    for (auto x : seq) util::write_gamma(builder, x);
    std::cout << "measured bits/int = " << (1.0 * builder.num_bits()) / seq.size() << std::endl;
    bit_vector_iterator bv_it(builder.data().data(), builder.data().size());
    for (auto x : seq) {
        uint64_t got = util::read_gamma(bv_it);
        REQUIRE_MESSAGE(got == x, "got " << got << " but expected " << x);
    }
}

TEST_CASE("delta") {
    const uint64_t u = rand() % 1000 + 1;  // at least 1
    std::cout << "u = " << u << std::endl;
    std::vector<uint64_t> seq = bits::test::get_sequence(sequence_length, u);
    bit_vector::builder builder;
    for (auto x : seq) util::write_delta(builder, x);
    std::cout << "measured bits/int = " << (1.0 * builder.num_bits()) / seq.size() << std::endl;
    bit_vector_iterator bv_it(builder.data().data(), builder.data().size());
    for (auto x : seq) {
        uint64_t got = util::read_delta(bv_it);
        REQUIRE_MESSAGE(got == x, "got " << got << " but expected " << x);
    }
}

TEST_CASE("rice") {
    const uint64_t u = rand() % 1000000 + 1;  // at least 1
    std::cout << "u = " << u << std::endl;
    std::vector<uint64_t> seq = bits::test::get_sequence(sequence_length, u);
    bit_vector::builder builder;
    for (uint64_t k = 1; k <= 32; ++k) {
        builder.clear();
        for (auto x : seq) util::write_rice(builder, x, k);
        std::cout << "measured bits/int = " << (1.0 * builder.num_bits()) / seq.size() << std::endl;
        bit_vector_iterator bv_it(builder.data().data(), builder.data().size());
        for (auto x : seq) {
            uint64_t got = util::read_rice(bv_it, k);
            REQUIRE_MESSAGE(got == x, "got " << got << " but expected " << x);
        }
    }
}
