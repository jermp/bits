#pragma once

#include "bit_vector.hpp"
#include "darray.hpp"
#include "compact_vector.hpp"

#include <iterator>

namespace bits {

template <bool index_zeros = false, bool encode_prefix_sum = false>
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
            n = n + 1;  // because I will add a zero at the beginning
        } else {
            if constexpr (std::is_same_v<typename std::iterator_traits<Iterator>::iterator_category,
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

        bit_vector_builder bvb_high_bits(n + (m_universe >> l) + 1);
        compact_vector::builder cv_builder_low_bits(n, l);

        uint64_t low_mask = (uint64_t(1) << l) - 1;
        uint64_t last = 0;

        // add a zero at the beginning
        if constexpr (encode_prefix_sum) {
            if (l) cv_builder_low_bits.push_back(0);
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
            if (l) cv_builder_low_bits.push_back(v & low_mask);
            bvb_high_bits.set((v >> l) + i, 1);
            last = v;
        }

        bit_vector(&bvb_high_bits).swap(m_high_bits);
        cv_builder_low_bits.build(m_low_bits);
        m_high_bits_d1.build(m_high_bits);
        if constexpr (index_zeros) m_high_bits_d0.build(m_high_bits);
    }

    struct iterator {
        iterator() : m_ef(nullptr) {}

        iterator(elias_fano const* ef, uint64_t pos = 0)
            : m_ef(ef), m_pos(pos), m_l(ef->m_low_bits.width()) {
            assert(m_pos < m_ef->size());
            assert(m_l < 64);
            if (m_ef->m_high_bits_d1.num_positions() == 0) return;
            uint64_t begin = m_ef->m_high_bits_d1.select(m_ef->m_high_bits, m_pos);
            m_high_enum = bit_vector::unary_iterator(m_ef->m_high_bits, begin);
            m_low_enum = m_ef->m_low_bits.at(m_pos);
        }

        bool good() const { return m_ef != nullptr; }
        bool has_next() const { return m_pos < m_ef->size(); }

        uint64_t next() {
            assert(good() and has_next());
            uint64_t high = m_high_enum.next();
            assert(high == m_ef->m_high_bits_d1.select(m_ef->m_high_bits, m_pos));
            uint64_t low = m_low_enum.value();
            uint64_t val = (((high - m_pos) << m_l) | low);
            ++m_pos;
            return val;
        }

    private:
        elias_fano const* m_ef;
        uint64_t m_pos;
        uint64_t m_l;
        bit_vector::unary_iterator m_high_enum;
        compact_vector::iterator m_low_enum;
    };

    iterator at(uint64_t pos) const {
        assert(pos < size());
        return iterator(this, pos);
    }

    iterator begin() const { return at(0); }
    iterator end() const { return at(size()); }

    inline uint64_t access(uint64_t i) const {
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
    inline uint64_t diff(uint64_t i) const {
        assert(i < size() && encode_prefix_sum);
        uint64_t low1 = m_low_bits.access(i);
        uint64_t low2 = m_low_bits.access(i + 1);
        uint64_t l = m_low_bits.width();
        uint64_t pos = m_high_bits_d1.select(m_high_bits, i);
        uint64_t h1 = pos - i;
        uint64_t h2 = bit_vector::unary_iterator(m_high_bits, pos + 1).next() - i - 1;
        uint64_t val1 = (h1 << l) | low1;
        uint64_t val2 = (h2 << l) | low2;
        return val2 - val1;
    }

    /*
        Return [position,value] of the leftmost element that is >= x.
        Return [size()-1,back()] if x >= back() (largest element).
    */
    inline std::pair<uint64_t, uint64_t> next_geq(uint64_t x) const { return next_geq_leftmost(x); }

    /*
        Return the position of the rightmost largest element <= x.
        Return size()-1 if x >= back() (largest element).
    */
    inline uint64_t prev_leq(uint64_t x) const {
        if (x >= back()) return size() - 1;
        auto [pos, val] = next_geq_rightmost(x);
        return pos - (val > x);
    }

    inline uint64_t back() const { return m_universe; }
    inline uint64_t size() const { return m_low_bits.size(); }

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
        Return [position,value] of the leftmost element that is >= x.
        Return [size()-1,back()] if x >= back() (largest element).
    */
    inline std::pair<uint64_t, uint64_t> next_geq_leftmost(uint64_t x) const {
        static_assert(index_zeros == true, "must build index on zeros");
        assert(m_high_bits_d0.num_positions());

        if (x >= back()) return {size() - 1, back()};

        uint64_t h_x = x >> m_low_bits.width();
        uint64_t begin = h_x ? m_high_bits_d0.select(m_high_bits, h_x - 1) - h_x + 1 : 0;
        assert(begin < size());

        // uint64_t end = m_high_bits_d0.select(m_high_bits, h_x) - h_x;
        // assert(end <= size());
        // assert(begin <= end);
        // return binary search for x in [begin, end)

        auto it = at(begin);
        uint64_t pos = begin;
        uint64_t val = it.next();
        while (val < x) {
            ++pos;
            val = it.next();
        }
        assert(val >= x);

        /* now pos is the position of the leftmost element that is >= x */
        assert(val == access(pos));
        return {pos, val};
    }

    /*
        Return [position,value] of the rightmost element that is >= x.
        Return [size(),back()] if x > back() (largest element).
    */
    inline std::pair<uint64_t, uint64_t> next_geq_rightmost(uint64_t x) const {
        auto [pos, val] = next_geq_leftmost(x);
        if (val == x and pos != size() - 1) {
            auto it = at(pos);
            val = it.next();
            while (val == x) {  // keep scanning to pick the rightmost one
                ++pos;
                val = it.next();
            }
            assert(val > x);
            return {pos - 1, x};
        }
        return {pos, val};
    }
};

}  // namespace bits