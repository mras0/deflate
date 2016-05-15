#ifndef DEFLATE_HUFFMAN_TABLE_H
#define DEFLATE_HUFFMAN_TABLE_H

#include "huffman_code.h"
#include <vector>

namespace deflate {

using huffman_table = std::vector<huffman_code>;

huffman_table make_huffman_table(const uint8_t* symbol_bit_lengths, const uint8_t* symbol_bit_lengths_end);
huffman_table make_huffman_table(const std::vector<uint8_t>& symbol_bit_lengths);
huffman_table make_default_huffman_table();
huffman_table make_default_huffman_len_table();

} // namespace deflate

#endif

