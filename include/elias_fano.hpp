#pragma once

#include "bit_vector.hpp"
#include "darray.hpp"
#include "compact_vector.hpp"

#include <iterator>

namespace bits {

template <  //
    /* build a succinct select index on the zeros of high_bits for efficient next_neq */
    bool index_zeros = false,  //
    /* if encode_prefix_sum = false, the sequence is assumed to be sorted */
    bool encode_prefix_sum = false  //
    >
struct elias_fano {
    elias_fano() : m_universe(0) {}

    template <typename Iterator>
    void encode(Iterator begin, uint64_t n) {
        if (n == 0) return;

        /* calculate the universe of representation, which is the largest element */
        assert(m_universe == 0);
        if constexpr (encode_prefix_sum) {
            auto tmp = begin;
            for (uint64_t i = 0; i != n; ++i, ++tmp) m_universe += *tmp;
            n = n + 1;  // because a zero is added at the beginning
        } else {
            if constexpr (std::is_same_v<typename Iterator::iterator_category,
                                         std::random_access_iterator_tag>) {
                m_universe = *(begin + (n - 1));
            } else {  // scan
                auto tmp = begin;
                for (uint64_t i = 0; i != n - 1; ++i, ++tmp)
                    ;
                m_universe = *tmp;
            }
        }

        /* This version takes at most: n*floor(log(U/n)) + 3*n bits */
        uint64_t l = uint64_t((n && m_universe / n) ? util::msb(m_universe / n) : 0);

        /* This version takes at most: n*ceil(log(U/n)) + 2*n bits */
        // uint64_t l = std::ceil(std::log2(static_cast<double>(m_universe) / n));

        /*
            Q. Which version is better?
            A. It depends on the indexes built on the high_bits.
        */

        bit_vector::builder bvb_high_bits(n + (m_universe >> l) + 1);
        compact_vector::builder cvb_low_bits(n, l);

        uint64_t low_mask = (uint64_t(1) << l) - 1;
        uint64_t last = 0;

        // add a zero at the beginning
        if constexpr (encode_prefix_sum) {
            if (l) cvb_low_bits.set(0, 0);
            bvb_high_bits.set(0, 1);
            n = n - 1;  // restore n
        }

        for (uint64_t i = 0; i != n; ++i, ++begin) {
            auto v = *begin;
            if constexpr (encode_prefix_sum) {
                v = v + last;             // prefix sum
            } else if (i and v < last) {  // check the order
                std::cerr << "error at " << i << "/" << n << ":\n";
                std::cerr << "last " << last << "\n";
                std::cerr << "current " << v << "\n";
                throw std::runtime_error("sequence is not sorted");
            }
            if (l) cvb_low_bits.set(i + encode_prefix_sum, v & low_mask);
            bvb_high_bits.set((v >> l) + i + encode_prefix_sum, 1);
            last = v;
        }

        bvb_high_bits.build(m_high_bits);
        cvb_low_bits.build(m_low_bits);
        m_high_bits_d1.build(m_high_bits);
        if constexpr (index_zeros) m_high_bits_d0.build(m_high_bits);
    }

    struct iterator {
        iterator() : m_ef(nullptr), m_pos(0), m_l(0), m_val(0) {}

        iterator(elias_fano const* ef, uint64_t pos = 0)
            : m_ef(ef)
            , m_pos(pos)
            , m_l(ef->m_low_bits.width())
            , m_val(0)  //
        {
            if (!has_next() or m_ef->m_high_bits_d1.num_positions() == 0) return;
            assert(m_l < 64);
            uint64_t begin = m_ef->m_high_bits_d1.select(m_ef->m_high_bits, m_pos);
            m_high_bits_it = m_ef->m_high_bits.get_iterator_at(begin);
            m_low_bits_it = m_ef->m_low_bits.get_iterator_at(m_pos);
            read_next_value();
        }

        bool has_next() const { return m_pos < m_ef->size(); }
        bool has_prev() const { return m_pos > 0; }
        uint64_t value() const { return m_val; }
        uint64_t position() const { return m_pos; }

        void next() {
            ++m_pos;
            if (!has_next()) return;
            read_next_value();
        }

        /*
            Return the value before the current position.
        */
        uint64_t prev_value() {
            assert(m_pos > 0);
            uint64_t pos = m_pos - 1;
            /*
                Read_next_value() sets the state ahead of 1 position,
                hence must go back by 2 to get previous value.
            */
            assert(m_high_bits_it.position() >= 2);
            uint64_t high = m_high_bits_it.prev(m_high_bits_it.position() - 2);
            assert(high == m_ef->m_high_bits_d1.select(m_ef->m_high_bits, pos));
            uint64_t low = *(m_low_bits_it - 2);
            return (((high - pos) << m_l) | low);
        }

    private:
        elias_fano const* m_ef;
        uint64_t m_pos;
        uint64_t m_l;
        uint64_t m_val;
        bit_vector::iterator m_high_bits_it;
        compact_vector::iterator m_low_bits_it;

        void read_next_value() {
            assert(m_pos < m_ef->size());
            uint64_t high = m_high_bits_it.next();
            assert(high == m_ef->m_high_bits_d1.select(m_ef->m_high_bits, m_pos));
            uint64_t low = *m_low_bits_it;
            m_val = (((high - m_pos) << m_l) | low);
            ++m_low_bits_it;
        }
    };

    iterator get_iterator_at(uint64_t pos) const { return iterator(this, pos); }
    iterator begin() const { return get_iterator_at(0); }

    uint64_t access(uint64_t i) const {
        assert(i < size());
        return ((m_high_bits_d1.select(m_high_bits, i) - i) << m_low_bits.width()) |
               m_low_bits.access(i);
    }

    /*
        If encode_prefix_sum = true, the iterator passed to the encode() method
        can yield non-sorted values, so that diff(i) returns the i-th element
        from the original sequence.
        Example. Assume the values are V = [3, 2, 5, 1, 16],
        then we will encode V' = [0, 3, 5, 10, 11, 27], so that
            diff(0) = V[0] = V'[1] - V'[0] = 3-0 = 3
            diff(1) = V[1] = V'[2] - V'[1] = 5-1 = 2
            diff(2) = V[2] = V'[3] - V'[2] = 10-5 = 5
            diff(3) = V[3] = V'[4] - V'[3] = 11-10 = 1
            diff(4) = V[4] = V'[5] - V'[4] = 27-11 = 16
    */
    uint64_t diff(uint64_t i) const {
        assert(i < size() && encode_prefix_sum);
        uint64_t low1 = m_low_bits.access(i);
        uint64_t low2 = m_low_bits.access(i + 1);
        uint64_t l = m_low_bits.width();
        uint64_t pos = m_high_bits_d1.select(m_high_bits, i);
        uint64_t h1 = pos - i;
        uint64_t h2 = m_high_bits.get_iterator_at(pos + 1).next() - i - 1;
        uint64_t val1 = (h1 << l) | low1;
        uint64_t val2 = (h2 << l) | low2;
        return val2 - val1;
    }

    /*
        Return [position,value] of the leftmost smallest element that is >= x.
        Return [size()-1,back()] if x > back().

        Example.

        1, 3, 3, 4, 5, 6, 6, 9, 12, 14, 17, 17
        0  1  2  3  4  5  6  7   8   9  10  11

        next_geq(0) = [0,1]
        next_geq(3) = [1,3]
        next_geq(6) = [5,6]
        next_geq(7) = [7,9]
        next_geq(17) = [10,17]
        next_geq(23) = [11,17] (saturate)
    */
    std::pair<uint64_t, uint64_t> next_geq(uint64_t x) const { return next_geq_leftmost(x); }

    /*
        Return the position of the rightmost largest element that is <= x.
        Return size()-1 if x >= back().
        Return -1 if x < front() (result is undefined).

        Example.

        1, 3, 3, 4, 5, 6, 6, 9, 12, 14, 17, 17
        0  1  2  3  4  5  6  7   8   9  10  11

        prev_leq(0) = -1 (undefined, because 0 < front() = 1)
        prev_leq(3) = 2
        prev_leq(6) = 6
        prev_leq(7) = 6
        prev_leq(17) = 11
        prev_leq(23) = 11 (saturate, because 23 >= back() = 17)
    */
    uint64_t prev_leq(uint64_t x) const {
        auto [pos, val] = next_geq_rightmost(x);
        return pos - (val > x);
    }

    uint64_t back() const { return m_universe; }
    uint64_t size() const { return m_low_bits.size(); }

    uint64_t num_bytes() const {
        return sizeof(m_universe) + m_high_bits.num_bytes() + m_high_bits_d1.num_bytes() +
               m_high_bits_d0.num_bytes() + m_low_bits.num_bytes();
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

private:
    uint64_t m_universe;
    bit_vector m_high_bits;
    darray1 m_high_bits_d1;
    darray0 m_high_bits_d0;
    compact_vector m_low_bits;

    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_universe);
        visitor.visit(t.m_high_bits);
        visitor.visit(t.m_high_bits_d1);
        visitor.visit(t.m_high_bits_d0);
        visitor.visit(t.m_low_bits);
    }

    /*
        Return [position,value] of the leftmost smallest element that is >= x.
        Return [size()-1,back()] if x > back().
    */
    std::pair<uint64_t, uint64_t> next_geq_leftmost(uint64_t x) const {
        static_assert(index_zeros == true, "must build index on zeros");
        assert(m_high_bits_d0.num_positions());

        if (x > back()) return {size() - 1, back()};

        uint64_t h_x = x >> m_low_bits.width();
        uint64_t begin = h_x ? m_high_bits_d0.select(m_high_bits, h_x - 1) - h_x + 1 : 0;
        assert(begin < size());

        // uint64_t end = m_high_bits_d0.select(m_high_bits, h_x) - h_x;
        // assert(end <= size());
        // assert(begin <= end);
        // return binary search for x in [begin, end)

        auto it = get_iterator_at(begin);
        uint64_t pos = begin;
        uint64_t val = it.value();
        while (val < x) {
            ++pos;
            /*
                Note: no need for bound checking here
                because x <= back(), hence pos cannot
                be equal to size().
            */
            it.next();
            val = it.value();
        }
        /* now pos is the position of the leftmost element that is >= x */
        assert(val >= x);
        assert(pos < size());
        assert(val == access(pos));
        return {pos, val};
    }

    /*
        Return [position,value] of the rightmost smallest element that is >= x.
        Return [size()-1,back()] if x >= back().
    */
    std::pair<uint64_t, uint64_t> next_geq_rightmost(uint64_t x) const {
        auto [pos, val] = next_geq_leftmost(x);
        if (val == x and pos != size() - 1) {
            auto it = get_iterator_at(pos);
            val = it.value();
            while (val == x) {  // keep scanning to pick the rightmost one
                ++pos;
                if (pos == size()) break;
                it.next();
                val = it.value();
            }
            assert(val >= x);
            return {pos - 1, x};
        }
        return {pos, val};
    }
};

}  // namespace bits