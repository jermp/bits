#pragma once

#include <cstddef>
#include <vector>

#include "essentials.hpp"
#include "util.hpp"

namespace bits {

struct bit_vector {
    struct builder {
        builder() { clear(); }

        builder(uint64_t num_bits) { resize(num_bits); }

        void clear() {
            m_num_bits = 0;
            m_data.clear();
            m_cur_word = nullptr;
        }

        void resize(uint64_t num_bits) {
            m_num_bits = num_bits;
            m_data.resize(essentials::words_for<uint64_t>(num_bits), 0);
        }

        void reserve(uint64_t num_bits) {
            m_data.reserve(essentials::words_for<uint64_t>(num_bits));
        }

        void build(bit_vector& bv) {
            bv.m_num_bits = m_num_bits;
            bv.m_data.swap(m_data);
            builder().swap(*this);
        }

        void swap(builder& other) {
            std::swap(m_num_bits, other.m_num_bits);
            std::swap(m_cur_word, other.m_cur_word);
            m_data.swap(other.m_data);
        }

        inline void push_back(bool b) {
            uint64_t pos_in_word = m_num_bits % 64;
            if (pos_in_word == 0) {
                m_data.push_back(0);
                m_cur_word = &m_data.back();
            }
            *m_cur_word |= (uint64_t)b << pos_in_word;
            ++m_num_bits;
        }

        inline void set(uint64_t pos, bool b = true) {
            assert(pos < num_bits());
            uint64_t word = pos >> 6;
            uint64_t pos_in_word = pos & 63;
            m_data[word] &= ~(uint64_t(1) << pos_in_word);
            m_data[word] |= uint64_t(b) << pos_in_word;
        }

        inline uint64_t get(uint64_t pos) const {
            assert(pos < num_bits());
            uint64_t word = pos >> 6;
            uint64_t pos_in_word = pos & 63;
            return m_data[word] >> pos_in_word & uint64_t(1);
        }

        inline void set_bits(uint64_t pos, uint64_t bits, uint64_t len) {
            assert(pos + len <= num_bits());
            // check there are no spurious bits
            assert(len == 64 || (bits >> len) == 0);
            if (!len) return;
            uint64_t mask = (len == 64) ? uint64_t(-1) : ((uint64_t(1) << len) - 1);
            uint64_t word = pos >> 6;
            uint64_t pos_in_word = pos & 63;

            m_data[word] &= ~(mask << pos_in_word);
            m_data[word] |= bits << pos_in_word;

            uint64_t stored = 64 - pos_in_word;
            if (stored < len) {
                m_data[word + 1] &= ~(mask >> stored);
                m_data[word + 1] |= bits >> stored;
            }
        }

        inline void append_bits(uint64_t bits, uint64_t len) {
            // check there are no spurious bits
            assert(len == 64 || (bits >> len) == 0);
            if (!len) return;
            uint64_t pos_in_word = m_num_bits & 63;
            m_num_bits += len;
            if (pos_in_word == 0) {
                m_data.push_back(bits);
            } else {
                *m_cur_word |= bits << pos_in_word;
                if (len > 64 - pos_in_word) m_data.push_back(bits >> (64 - pos_in_word));
            }
            m_cur_word = &m_data.back();
        }

        inline uint64_t get_word64(uint64_t pos) const {
            assert(pos < num_bits());
            uint64_t block = pos >> 6;
            uint64_t shift = pos & 63;
            uint64_t word = m_data[block] >> shift;
            if (shift && block + 1 < m_data.size()) word |= m_data[block + 1] << (64 - shift);
            return word;
        }

        void append(builder const& rhs) {
            if (!rhs.num_bits()) return;

            uint64_t pos = m_data.size();
            uint64_t shift = num_bits() % 64;
            m_num_bits = num_bits() + rhs.num_bits();
            m_data.resize(essentials::words_for<uint64_t>(m_num_bits));

            if (shift == 0) {  // word-aligned, easy case
                std::copy(rhs.m_data.begin(), rhs.m_data.end(), m_data.begin() + ptrdiff_t(pos));
            } else {
                uint64_t* cur_word = &m_data.front() + pos - 1;
                for (uint64_t i = 0; i < rhs.m_data.size() - 1; ++i) {
                    uint64_t w = rhs.m_data[i];
                    *cur_word |= w << shift;
                    *++cur_word = w >> (64 - shift);
                }
                *cur_word |= rhs.m_data.back() << shift;
                if (cur_word < &m_data.back()) *++cur_word = rhs.m_data.back() >> (64 - shift);
            }
            m_cur_word = &m_data.back();
        }

        uint64_t num_bits() const { return m_num_bits; }
        std::vector<uint64_t> const& data() const { return m_data; }

    private:
        uint64_t m_num_bits;
        uint64_t* m_cur_word;
        std::vector<uint64_t> m_data;
    };

    bit_vector() : m_num_bits(0) {}

    void swap(bit_vector& other) {
        std::swap(other.m_num_bits, m_num_bits);
        other.m_data.swap(m_data);
    }

    // get i-th bit
    inline uint64_t operator[](uint64_t i) const {
        assert(i < num_bits());
        uint64_t block = i >> 6;
        uint64_t shift = i & 63;
        return m_data[block] >> shift & uint64_t(1);
    }

    inline uint64_t get_bits(uint64_t pos, uint64_t len) const {
        assert(pos + len <= num_bits());
        if (!len) return 0;
        uint64_t block = pos >> 6;
        uint64_t shift = pos & 63;
        uint64_t mask = -(len == 64) | ((1ULL << len) - 1);
        if (shift + len <= 64) {
            return m_data[block] >> shift & mask;
        } else {
            return (m_data[block] >> shift) | (m_data[block + 1] << (64 - shift) & mask);
        }
    }

    // fast and unsafe version: it retrieves at least 56 bits
    inline uint64_t get_word56(uint64_t pos) const {
        const char* base_ptr = reinterpret_cast<const char*>(m_data.data());
        return *(reinterpret_cast<uint64_t const*>(base_ptr + (pos >> 3))) >> (pos & 7);
    }

    // pad with zeros if extension further size is needed
    inline uint64_t get_word64(uint64_t pos) const {
        assert(pos < num_bits());
        uint64_t block = pos >> 6;
        uint64_t shift = pos & 63;
        uint64_t word = m_data[block] >> shift;
        if (shift && block + 1 < m_data.size()) { word |= m_data[block + 1] << (64 - shift); }
        return word;
    }

    struct unary_iterator {
        unary_iterator() : m_data(0), m_position(0), m_buf(0) {}

        unary_iterator(bit_vector const& bv, uint64_t pos = 0) {
            m_data = bv.data().data();
            m_position = pos;
            m_buf = m_data[pos >> 6];
            // clear low bits
            m_buf &= uint64_t(-1) << (pos & 63);
        }

        uint64_t position() const { return m_position; }

        uint64_t next() {
            uint64_t pos_in_word;
            uint64_t buf = m_buf;
            while (!util::lsbll(buf, pos_in_word)) {
                m_position += 64;
                buf = m_data[m_position >> 6];
            }

            m_buf = buf & (buf - 1);  // clear LSB
            m_position = (m_position & ~uint64_t(63)) + pos_in_word;
            return m_position;
        }

        // skip to the k-th one after the current position
        void skip(uint64_t k) {
            uint64_t skipped = 0;
            uint64_t buf = m_buf;
            uint64_t w = 0;
            while (skipped + (w = util::popcount(buf)) <= k) {
                skipped += w;
                m_position += 64;
                buf = m_data[m_position / 64];
            }
            assert(buf);
            uint64_t pos_in_word = util::select_in_word(buf, k - skipped);
            m_buf = buf & (uint64_t(-1) << pos_in_word);
            m_position = (m_position & ~uint64_t(63)) + pos_in_word;
        }

        // skip to the k-th zero after the current position
        void skip0(uint64_t k) {
            uint64_t skipped = 0;
            uint64_t pos_in_word = m_position % 64;
            uint64_t buf = ~m_buf & (uint64_t(-1) << pos_in_word);
            uint64_t w = 0;
            while (skipped + (w = util::popcount(buf)) <= k) {
                skipped += w;
                m_position += 64;
                buf = ~m_data[m_position / 64];
            }
            assert(buf);
            pos_in_word = util::select_in_word(buf, k - skipped);
            m_buf = ~buf & (uint64_t(-1) << pos_in_word);
            m_position = (m_position & ~uint64_t(63)) + pos_in_word;
        }

    private:
        uint64_t const* m_data;
        uint64_t m_position;
        uint64_t m_buf;
    };

    uint64_t num_bits() const { return m_num_bits; }
    std::vector<uint64_t> const& data() const { return m_data; }

    uint64_t num_bytes() const { return sizeof(m_num_bits) + essentials::vec_bytes(m_data); }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

protected:
    uint64_t m_num_bits;
    std::vector<uint64_t> m_data;

    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_num_bits);
        visitor.visit(t.m_data);
    }
};

struct bit_vector_iterator {
    bit_vector_iterator() : m_data(nullptr), m_num_64bit_words(0), m_pos(0), m_buf(0), m_avail(0) {}

    bit_vector_iterator(uint64_t const* data, uint64_t num_64bit_words, uint64_t pos = 0)
        : m_data(data), m_num_64bit_words(num_64bit_words) {
        at(pos);
    }

    void at(uint64_t pos) {
        m_pos = pos;
        m_buf = 0;
        m_avail = 0;
    }

    void at_and_clear_low_bits(uint64_t pos) {
        m_pos = pos;
        m_buf = m_data[pos / 64];
        m_buf &= uint64_t(-1) << (m_pos & 63);  // clear low bits
    }

    /*
        Return the next l bits from the current position and advance by l bits.
    */
    inline uint64_t take(uint64_t l) {
        assert(l <= 64);
        if (m_avail < l) fill_buf();
        uint64_t val;
        if (l != 64) {
            val = m_buf & ((uint64_t(1) << l) - 1);
            m_buf >>= l;
        } else {
            val = m_buf;
        }
        m_avail -= l;
        m_pos += l;
        return val;
    }

    uint64_t next() {
        uint64_t pos_in_word;
        uint64_t buf = m_buf;
        while (!util::lsbll(buf, pos_in_word)) {
            m_pos += 64;
            buf = m_data[m_pos >> 6];
        }
        m_buf = buf & (buf - 1);  // clear LSB
        m_pos = (m_pos & ~uint64_t(63)) + pos_in_word;
        return m_pos;
    }

    /*
        Skip all zeros from the current position and
        return the number of skipped zeros.
    */
    inline uint64_t skip_zeros() {
        uint64_t zeros = 0;
        while (m_buf == 0) {
            m_pos += m_avail;
            zeros += m_avail;
            fill_buf();
        }
        uint64_t l = util::lsbll(m_buf);
        m_buf >>= l;
        m_buf >>= 1;
        m_avail -= l + 1;
        m_pos += l + 1;
        return zeros + l;
    }

    inline uint64_t position() const { return m_pos; }

    inline void fill_buf() {
        m_buf = get_next_word64();
        m_avail = 64;
    }

private:
    uint64_t get_next_word64() {
        uint64_t block = m_pos / 64;
        uint64_t shift = m_pos % 64;
        uint64_t word = m_data[block] >> shift;
        if (shift && block + 1 < m_num_64bit_words) { word |= m_data[block + 1] << (64 - shift); }
        return word;
    }

    uint64_t const* m_data;
    uint64_t m_num_64bit_words;
    uint64_t m_pos;
    uint64_t m_buf;
    uint64_t m_avail;
};

}  // namespace bits