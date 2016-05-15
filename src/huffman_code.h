#ifndef DEFLATE_HUFFMAN_CODE_H
#define DEFLATE_HUFFMAN_CODE_H

#include <stdint.h>

namespace deflate {

static constexpr int max_bits = 15;

struct huffman_code {
    uint8_t  len;
    uint32_t value;

    bool valid() const {
        return len > 0 && len <= max_bits && len <= 8 * sizeof(value) && (value >> len) == 0;
    }
};

inline bool operator==(const huffman_code& l, const huffman_code& r) { return l.len == r.len && l.value == r.value; }
inline bool operator!=(const huffman_code& l, const huffman_code& r) { return !(l == r); }

} // namespace deflate

#endif