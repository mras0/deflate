#include "huffman_table.h"
#include <algorithm>
#include <cassert>

namespace deflate {

huffman_table make_huffman_table(const uint8_t* symbol_bit_lengths, const uint8_t* symbol_bit_lengths_end)
{
    const auto max_bit_length = *std::max_element(symbol_bit_lengths, symbol_bit_lengths_end);
    assert(max_bit_length <= max_bits);
    // Count the number of codes for each code length. Let bl_count[N] be the number of codes of length N, N >= 1.
    std::vector<int> bl_count(max_bit_length+1);
    for (auto bl = symbol_bit_lengths; bl != symbol_bit_lengths_end; ++bl) {
        if (*bl > 0) {
            ++bl_count[*bl];
        }
    }

    //  Find the numerical value of the smallest code for each code length
    uint16_t code = 0;
    assert(bl_count[0] == 0);
    std::vector<uint16_t> next_code(max_bit_length+1);
    for (int bits = 1; bits <= max_bit_length; bits++) {
        code = static_cast<uint16_t>((code + bl_count[bits-1]) << 1);
        next_code[bits] = code;
    }

    // Assign numerical values to all codes, using consecutive values for all codes of the same length with the base
    // values determined at step 2. Codes that are never used (which have a bit length of zero) must not be assigned
    // a value.

    huffman_table codes(symbol_bit_lengths_end - symbol_bit_lengths);
    for (int n = 0; n < (int)codes.size(); n++) {
        const auto len = symbol_bit_lengths[n];
        if (len != 0) {            
            codes[n] = huffman_code{len, next_code[len]};
            next_code[len]++;
        }
    }
    return codes;
}

huffman_table make_huffman_table(const std::vector<uint8_t>& symbol_bit_lengths)
{
    return make_huffman_table(symbol_bit_lengths.data(), symbol_bit_lengths.data() + symbol_bit_lengths.size());
}

huffman_table make_default_huffman_table()
{
    /*
    Lit Value    Bits        Codes
    ---------    ----        -----
    0 - 143     8          00110000 through
    10111111
    144 - 255     9          110010000 through
    111111111
    256 - 279     7          0000000 through
    0010111
    280 - 287     8          11000000 through
    11000111
    */
    constexpr int num_symbols = 288;
    std::vector<uint8_t> sbl(num_symbols);
    int i = 0;
    while (i < 144) sbl[i++] = 8;
    while (i < 256) sbl[i++] = 9;
    while (i < 280) sbl[i++] = 7;
    while (i < num_symbols) sbl[i++] = 8;
    return make_huffman_table(sbl);
}

huffman_table make_default_huffman_len_table()
{
    // Distance codes 0-31 are represented by (fixed-length) 5-bit codes
    std::vector<uint8_t> length_bit_lengths(32, 5);
    return make_huffman_table(length_bit_lengths);
}

} // namespace deflate
