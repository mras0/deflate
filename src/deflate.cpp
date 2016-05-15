#include "deflate.h"
#include "huffman_tree.h"
#include "huffman_table.h"

#include <stdexcept>

namespace deflate {

enum class block_type { uncompressed, fixed_huffman, dynamic_huffman, reserved };

int decode(const huffman_tree& t, bit_stream& bs)
{
    int value = huffman_tree::max_symbols;
    if (bs.potentially_available_bits() >= t.table_bits()) {
        bs.ensure_bits(t.table_bits());
        auto te = t.next_from_bits(bs.peek_bits(t.table_bits()), t.table_bits());
        bs.consume_bits(te.len);
        value = te.index;
    }
    while (value >= huffman_tree::max_symbols) {
        value = t.branch(value - huffman_tree::max_symbols, !!bs.get_bit());
    }
    return value;
}

std::vector<uint8_t> deflate(bit_stream& bs)
{
    enum alphabet {
        lit_min      = 0,
        lit_max      = 255,
        end_of_block = 256,
        len_min      = 257,
        len_max      = 285,
    };
    constexpr int extra_bits[1+len_max - len_min] = {
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
    };
    constexpr int lengths[1+len_max-len_min] = {
        3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258,
    };

    constexpr int distance_extra_bits[32] = {
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13,
    };
    constexpr int distance_length[32] = {
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
    };

    auto invalid = [] () { assert(false); throw std::runtime_error("Invalid deflate stream"); };

    static const auto default_lit_len_tree = make_huffman_tree(make_default_huffman_table(), 9);
    static const auto default_dist_tree    = make_huffman_tree(make_default_huffman_len_table(), 5);

    std::vector<uint8_t> output;
    bool last_block = false;
    while (!last_block) {
        // Block header bits
        last_block = !!bs.get_bit();
        const auto type  = static_cast<block_type>(bs.get_bits(2));
        //std::cout << "final " << last_block << " type " << type << std::endl;
        assert(type != block_type::reserved);

        if (type == block_type::uncompressed) {
            assert(false);
            // skip any remaining bits in current partially
            //     processed byte
            //     read LEN and NLEN (each 16-bits)
            //     copy LEN bytes of data to output
        } else {
            huffman_tree lit_len_tree;
            huffman_tree dist_tree;

            if (type == block_type::dynamic_huffman) {
                // read representation of code trees
                const int hlit  = len_min + bs.get_bits(5); // 5 Bits: HLIT, # of Literal/Length codes - 257 (257 - 286)
                const int hdist = 1 + bs.get_bits(5);       // 5 Bits: HDIST, # of Distance codes - 1        (1 - 32)
                const int hclen = 4 + bs.get_bits(4);       // 4 Bits: HCLEN, # of Code Length codes - 4     (4 - 19)

                                                            //std::cout << "hlit  = " << hlit  << "\n";
                                                            //std::cout << "hdist = " << hdist << "\n";
                                                            //std::cout << "hclen = " << hclen << "\n";

                constexpr int max_code_lengths = 19;
                uint8_t code_lengths[max_code_lengths]={0};
                for (int i = 0; i < hclen; ++i) {
                    constexpr int alphabet_permute[max_code_lengths] = {
                        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
                    };
                    code_lengths[alphabet_permute[i]] = static_cast<uint8_t>(bs.get_bits(3));
                }

                huffman_tree cl_tree = make_huffman_tree(make_huffman_table(code_lengths, code_lengths+max_code_lengths), 7);

                std::vector<uint8_t> cl2(hlit + hdist);
                for (int i = 0; i < hlit + hdist;) {
                    auto cl_val = static_cast<uint8_t>(decode(cl_tree, bs));
                    int count = 1;
                    if (cl_val <= 15) {
                        // 0 - 15: Represent code lengths of 0 - 15
                    } else if (cl_val == 16) {
                        // 16: Copy the previous code length 3 - 6 times.
                        if (i == 0) invalid();
                        cl_val = cl2[i-1];
                        // The next 2 bits indicate repeat length
                        // (0 = 3, ... , 3 = 6)
                        // Example:  Codes 8, 16 (+2 bits 11),
                        // 16 (+2 bits 10) will expand to
                        // 12 code lengths of 8 (1 + 6 + 5)
                        count = 3 + bs.get_bits(2);
                    } else if (cl_val == 17) {
                        // 17: Repeat a code length of 0 for 3 - 10 times.
                        // (3 bits of length)
                        cl_val = 0;
                        count = 3 + bs.get_bits(3);
                    } else {
                        // 18: Repeat a code length of 0 for 11 - 138 times
                        // (7 bits of length)
                        assert(cl_val == 18);
                        cl_val = 0;
                        count = 11 + bs.get_bits(7);
                    }
                    assert(cl_val <= max_bits);
                    if (static_cast<size_t>(count + i) > cl2.size()) invalid();
                    while (count--) {
                        //std::cout << (i > hlit ? "dist" : "lit") <<  "[" << (i > hlit ? i - hlit : i) << "] = " << (int)cl_val << std::endl;
                        cl2[i++] = cl_val;
                    }
                }

                lit_len_tree = make_huffman_tree(make_huffman_table(cl2.data(), cl2.data() + hlit), 9);
                dist_tree    = make_huffman_tree(make_huffman_table(cl2.data() + hlit, cl2.data() + hlit + hdist), 6);

                //static int n = 0;
                //write_huff_tree("lit"+std::to_string(n)+".gv", lit_len_tree);
                //write_huff_tree("dist"+std::to_string(n)+".gv", dist_tree);
                //++n;
            } else {
                assert(type == block_type::fixed_huffman);
                lit_len_tree = default_lit_len_tree;
                dist_tree    = default_dist_tree;
            }


            for (;;) {
                // decode literal/length value from input stream
                const int value = decode(lit_len_tree, bs);
                if (value <= lit_max) {
                    // if value < 256
                    //    copy value (literal byte) to output stream
                    output.push_back(static_cast<uint8_t>(value));
                    //std::cout << "Lit " << litrep(value) << "\n";
                } else if (value == end_of_block) {
                    // if value = end of block (256)
                    //    break from loop
                    break;
                } else {
                    // otherwise (value = 257..285)
                    //    decode distance from input stream
                    //
                    //    move backwards distance bytes in the output
                    //    stream, and copy length bytes from this
                    //    position to the output stream.
                    const auto eb = extra_bits[value - len_min];
                    int len = lengths[value - len_min];
                    if (eb) {
                        len += bs.get_bits(eb);
                    }

                    int dist = decode(dist_tree, bs);
                    const int dist_extra_bits = distance_extra_bits[dist];
                    int dist_bytes = distance_length[dist];
                    if (dist_extra_bits) {
                        dist_bytes += bs.get_bits(dist_extra_bits);
                    }

                    //std::cout << "<" << len << ", " << dist_bytes << ">" << std::endl;
                    if (static_cast<size_t>(dist_bytes) > output.size()) invalid();
                    const int old_size = static_cast<int>(output.size());
                    const int src = old_size - dist_bytes;
                    output.resize(old_size + len);
                    for (int i = 0; i < len; ++i) {
                        output[old_size+i] = output[src+i];
                    }
                }
            }
        }
    }
    return output;
}

} // namespace deflate
