#pragma once

#include "util.hpp"
#include "bit_vector.hpp"

namespace bits {

/*
    The following class implements an index on top of an uncompressed
    bitvector B[0..n) to support select queries.

    The solution implemented here is described in

        Okanohara, Daisuke, and Kunihiko Sadakane. 2007.
        Practical entropy-compressed rank/select dictionary.
        In Proceedings of the 9-th Workshop on Algorithm Engineering
        and Experiments (ALENEX), pages 60-70.

    under the name "darray" (short for "dense" array).

    The bitvector is split into variable-length blocks, each
    containing L ones (except for, possibly, the last block).
    A block is said to be "sparse" if its length is >= L2;
    otherwise it is said "dense". Sparse blocks are represented
    verbatim, i.e., the positions of the L ones are coded using 64-bit
    integers. A dense block, instead, is sparsified: we keep
    one position every L3 positions. The positions are coded relatively
    to the beginning of each block, hence using log2(L2) bits per
    position.

    Let m be the number of bits set in B. The data structure stores:
    - An array `block_inventory[0..m/L)` such that `block_inventory[i]`
      holds the result of Select(iL) if block i is dense; it holds the start position
      in `overflow_positions` of the 1-bit positions of the block i if it is sparse.
      Its space is m/L*64 bits.
    - An array `overflow_positions` holding the positions of the L ones
      in sparse blocks. As we have at most m/L2 sparse blocks, its space is
      m/L2*L*64 bits at most.
    - An array `subblock_inventory[0..m/L3)` such that `subblock_inventory[i]`
      holds the result of Select(iL3). Its space is m/L3*log2(L2) bits.

    A Select(i) query first checks if the block i/L is sparse by accessing
    `block_inventory[i/L]`: if it is, then the query is solved in O(1) as
    `overflow_positions[block_inventory[i/L] + i%L]`; otherwise the corresponding
    sub-block is accessed and a sequential scan of at most L2 bits is performed.
    If p is the position computed during the scan, the final result is
    `block_inventory[i/L] + subblock_inventory[i/L3] + p`.

    If Scan(Q) = Q/B is the number of cache-misses involved in a sequential scan
    of Q bits with a cache-line of B bits, it follows that the number of cache-misses
    per Select(i) query are:
        2, if position i belongs to a sparse block;
        2 + Scan(L2), if i belongs to a dense block.

    This implementation uses, by default:

    L =  1,024 (block_size)
    L2 = 65,536 (so that each position in a dense block can be coded
                 using 16-bit integers)
    L3 = 32 (subblock_size)

    For these block sizes, we have a space usage of at most:
    m/2^10*64 (block_inventory) +
    m/2^16*2^10*64 + (sparse blocks) +
    m/2^5*2^4 (dense blocks) =
    25/16 m = 1.5625 m bits.

    (When used to index the high bitvector of Elias-Fano, sparse blocks are rare;
     so the space usage is likely close to 9/16 m.)

    If 0.0 <= d = m/n <= 1.0 is the density of the bitvector, the index costs at most
    25/16 dn extra bits, for a total of n(1+25/16d) bits.
*/

template <                       //
    typename WordGetter,         //
    uint64_t block_size = 1024,  //
    uint64_t subblock_size = 32  //
    >
struct darray {
    darray() : m_positions(0) {}

    void build(bit_vector const& B) {
        std::vector<uint64_t> const& data = B.data();
        std::vector<uint64_t> cur_block_positions;
        std::vector<int64_t> block_inventory;
        std::vector<uint16_t> subblock_inventory;
        std::vector<uint64_t> overflow_positions;

        for (uint64_t word_idx = 0; word_idx < data.size(); ++word_idx) {
            uint64_t cur_pos = word_idx << 6;
            uint64_t cur_word = WordGetter()(data, word_idx);
            uint64_t l;
            while (util::lsbll(cur_word, l)) {
                cur_pos += l;
                cur_word >>= l;
                if (cur_pos >= B.num_bits()) break;

                cur_block_positions.push_back(cur_pos);

                if (cur_block_positions.size() == block_size) {
                    flush_cur_block(cur_block_positions, block_inventory, subblock_inventory,
                                    overflow_positions);
                }

                // can't do >>= l + 1, can be 64
                cur_word >>= 1;
                cur_pos += 1;
                m_positions += 1;
            }
        }
        if (cur_block_positions.size()) {
            flush_cur_block(cur_block_positions, block_inventory, subblock_inventory,
                            overflow_positions);
        }
        m_block_inventory.swap(block_inventory);
        m_subblock_inventory.swap(subblock_inventory);
        m_overflow_positions.swap(overflow_positions);

        // {
        //     std::cout << "num_blocks = " << m_block_inventory.size()
        //               << " (expected = " << ((m_positions + block_size - 1) / block_size) << ")"
        //               << std::endl;
        //     std::cout << "num_subblocks = " << m_subblock_inventory.size()
        //               << " (expected = " << ((m_positions + subblock_size - 1) / subblock_size)
        //               << ")" << std::endl;

        //     std::cout << "block_inventory: got "
        //               << (essentials::vec_bytes(m_block_inventory) * 8.0 - 64) / m_positions
        //               << " bits/bit (expected "
        //               << ((m_positions + block_size - 1) / block_size * 64.0) / m_positions <<
        //               ")"
        //               << std::endl;

        //     std::cout << "overflow_positions: got "
        //               << (essentials::vec_bytes(m_overflow_positions) * 8.0 - 64) / m_positions
        //               << " bits/bit (at most 1)" << std::endl;

        //     std::cout << "subblock_inventory: got "
        //               << (essentials::vec_bytes(m_subblock_inventory) * 8.0 - 16) / m_positions
        //               << " bits/bit (expected "
        //               << ((m_positions + subblock_size - 1) / subblock_size * 16.0) / m_positions
        //               << ")" << std::endl;
        // }
    }

    /*
        Return the position of the i-th bit set in B,
        for any 0 <= i < num_positions();
    */
    inline uint64_t select(bit_vector const& B, uint64_t i) const {
        assert(i < num_positions());
        uint64_t block = i / block_size;
        int64_t block_pos = m_block_inventory[block];
        if (block_pos < 0) {  // sparse block
            uint64_t overflow_pos = uint64_t(-block_pos - 1);
            return m_overflow_positions[overflow_pos + (i & (block_size - 1))];
        }

        uint64_t subblock = i / subblock_size;
        uint64_t start_pos = uint64_t(block_pos) + m_subblock_inventory[subblock];
        uint64_t reminder = i & (subblock_size - 1);
        if (!reminder) return start_pos;

        std::vector<uint64_t> const& data = B.data();
        uint64_t word_idx = start_pos >> 6;
        uint64_t word_shift = start_pos & 63;
        uint64_t word = WordGetter()(data, word_idx) & (uint64_t(-1) << word_shift);
        while (true) {
            uint64_t popcnt = util::popcount(word);
            if (reminder < popcnt) break;
            reminder -= popcnt;
            word = WordGetter()(data, ++word_idx);
        }
        return (word_idx << 6) + util::select_in_word(word, reminder);
    }

    inline uint64_t num_positions() const { return m_positions; }

    uint64_t num_bytes() const {
        return sizeof(m_positions) + essentials::vec_bytes(m_block_inventory) +
               essentials::vec_bytes(m_subblock_inventory) +
               essentials::vec_bytes(m_overflow_positions);
    }

    void swap(darray& other) {
        std::swap(m_positions, other.m_positions);
        m_block_inventory.swap(other.m_block_inventory);
        m_subblock_inventory.swap(other.m_subblock_inventory);
        m_overflow_positions.swap(other.m_overflow_positions);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

protected:
    uint64_t m_positions;
    std::vector<int64_t> m_block_inventory;
    std::vector<uint16_t> m_subblock_inventory;
    std::vector<uint64_t> m_overflow_positions;

    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_positions);
        visitor.visit(t.m_block_inventory);
        visitor.visit(t.m_subblock_inventory);
        visitor.visit(t.m_overflow_positions);
    }

    static void flush_cur_block(std::vector<uint64_t>& cur_block_positions,
                                std::vector<int64_t>& block_inventory,
                                std::vector<uint16_t>& subblock_inventory,
                                std::vector<uint64_t>& overflow_positions)  //
    {
        if (cur_block_positions.back() - cur_block_positions.front() < (1ULL << 16))  // dense case
        {
            block_inventory.push_back(int64_t(cur_block_positions.front()));
            for (uint64_t i = 0; i < cur_block_positions.size(); i += subblock_size) {
                subblock_inventory.push_back(
                    uint16_t(cur_block_positions[i] - cur_block_positions.front()));
            }
        } else  // sparse case
        {
            block_inventory.push_back(-int64_t(overflow_positions.size()) - 1);
            for (uint64_t i = 0; i < cur_block_positions.size(); ++i) {
                overflow_positions.push_back(cur_block_positions[i]);
            }
            for (uint64_t i = 0; i < cur_block_positions.size(); i += subblock_size) {
                subblock_inventory.push_back(uint16_t(-1));
            }
        }
        cur_block_positions.clear();
    }
};

namespace util {

struct identity_getter {
    uint64_t operator()(std::vector<uint64_t> const& data, uint64_t i) const { return data[i]; }
};

struct negating_getter {
    uint64_t operator()(std::vector<uint64_t> const& data, uint64_t i) const { return ~data[i]; }
};

}  // namespace util

typedef darray<util::identity_getter> darray1;  // take positions of 1s
typedef darray<util::negating_getter> darray0;  // take positions of 0s

}  // namespace bits