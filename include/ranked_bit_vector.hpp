#pragma once

#include "bit_vector.hpp"

namespace bits {

struct ranked_bit_vector : public bit_vector {
    ranked_bit_vector() : bit_vector() {}

    void build_index() {
        std::vector<uint64_t> block_rank_pairs;
        uint64_t next_rank = 0;
        uint64_t cur_subrank = 0;
        uint64_t subranks = 0;
        block_rank_pairs.push_back(0);
        for (uint64_t i = 0; i < m_data.size(); ++i) {
            uint64_t word_pop = util::popcount(m_data[i]);
            uint64_t shift = i % block_size;
            if (shift) {
                subranks <<= 9;
                subranks |= cur_subrank;
            }
            next_rank += word_pop;
            cur_subrank += word_pop;

            if (shift == block_size - 1) {
                block_rank_pairs.push_back(subranks);
                block_rank_pairs.push_back(next_rank);
                subranks = 0;
                cur_subrank = 0;
            }
        }
        uint64_t left = block_size - m_data.size() % block_size;
        for (uint64_t i = 0; i < left; ++i) {
            subranks <<= 9;
            subranks |= cur_subrank;
        }
        block_rank_pairs.push_back(subranks);

        if (m_data.size() % block_size) {
            block_rank_pairs.push_back(next_rank);
            block_rank_pairs.push_back(0);
        }

        m_block_rank_pairs.swap(block_rank_pairs);
    }

    inline uint64_t num_ones() const { return *(m_block_rank_pairs.end() - 2); }
    inline uint64_t num_zeros() const { return num_bits() - num_ones(); }

    /*
        Return the number of ones in B[0..i).
    */
    inline uint64_t rank1(uint64_t i) const {
        assert(i <= num_bits());
        if (i == num_bits()) return num_ones();
        uint64_t sub_block = i / 64;
        uint64_t r = sub_block_rank(sub_block);
        uint64_t sub_left = i % 64;
        if (sub_left) r += util::popcount(m_data[sub_block] << (64 - sub_left));
        return r;
    }

    /*
        Return the number of zeros in B[0..i).
    */
    inline uint64_t rank0(uint64_t i) const {
        assert(i <= num_bits());
        return i - rank1(i);
    }

    uint64_t num_bytes() const {
        return bit_vector::num_bytes() + essentials::vec_bytes(m_block_rank_pairs);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        bit_vector::visit(visitor);
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        bit_vector::visit(visitor);
        visit_impl(visitor, *this);
    }

protected:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_block_rank_pairs);
    }

    inline uint64_t block_rank(uint64_t block) const { return m_block_rank_pairs[block * 2]; }

    inline uint64_t sub_block_rank(uint64_t sub_block) const {
        uint64_t r = 0;
        uint64_t block = sub_block / block_size;
        r += block_rank(block);
        uint64_t left = sub_block % block_size;
        r += sub_block_ranks(block) >> ((7 - left) * 9) & 0x1FF;
        return r;
    }

    inline uint64_t sub_block_ranks(uint64_t block) const {
        return m_block_rank_pairs[block * 2 + 1];
    }

    static const uint64_t block_size = 8;  // in 64bit words
    std::vector<uint64_t> m_block_rank_pairs;
};

}  // namespace bits