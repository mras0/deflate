#ifndef DEFLATE_BIT_STREAM_H
#define DEFLATE_BIT_STREAM_H

#include <stdint.h>
#include <cassert>

namespace deflate {
/*
* Data elements are packed into bytes in order of
increasing bit number within the byte, i.e., starting
with the least-significant bit of the byte.
* Data elements other than Huffman codes are packed
starting with the least-significant bit of the data
element.
* Huffman codes are packed starting with the most-
significant bit of the code.
*/

class bit_stream {
public:
    explicit bit_stream(const uint8_t* begin, const uint8_t* end) : data_(begin), len_(static_cast<int>(end - begin)) {
    }

    template<int size>
    explicit bit_stream(const uint8_t (&data)[size]) : data_(data), len_(size) {
    }

    void ensure_bits(int num_bits) {
        assert(num_bits > 0 && num_bits <= 16);
        while (avail_ < num_bits) {
            assert(pos_ < len_);
            const auto b = data_[pos_++];
            bits_   = (static_cast<uint32_t>(b) << avail_) | bits_;
            avail_ += 8;
        }
    }

    int potentially_available_bits() const {
        const auto remaining = len_ - pos_;
        return remaining >= 2 ? 16 : avail_ + 8 * remaining;
    }

    uint32_t peek_bits(int num_bits) {
        assert(num_bits > 0 && num_bits <= avail_);
        return bits_ & ((1 << num_bits) - 1);
    }

    void consume_bits(int num_bits) {
        assert(num_bits > 0 && num_bits <= avail_);
        avail_ -= num_bits;
        bits_ >>= num_bits;
    }

    int available_bits() const {
        return avail_;
    }

    uint8_t get_bit() {
        return static_cast<uint8_t>(get_bits(1));
    }

    uint32_t get_bits(int num_bits) {
        ensure_bits(num_bits);
        const auto res = peek_bits(num_bits);
        consume_bits(num_bits);
        return res;
    }

private:
    const uint8_t*  data_;
    int             len_;
    int             pos_   = 0;
    uint32_t        bits_  = 0;
    int             avail_ = 0;
};

} // namespace deflate

#endif
