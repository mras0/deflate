#include "deflate.h"
#include "huffman_tree.h"
#include "huffman_table.h"

#include <stdexcept>
#include <string.h>

namespace deflate {

enum class block_type { uncompressed, fixed_huffman, dynamic_huffman, reserved };

constexpr int max_match_length = 258;

int decode(const huffman_tree& t, bit_stream& bs)
{
    int value = huffman_tree::max_symbols;
    if (bs.potentially_available_bits() >= t.table_bits()) {
        bs.ensure_bits(t.table_bits());
        auto te = t.next_from_bits(bs.peek_bits(t.table_bits()), t.table_bits());
        bs.consume_bits(te.len());
        value = te.index();
    }
    while (value >= huffman_tree::max_symbols) {
        value = t.branch(value - huffman_tree::max_symbols, !!bs.get_bit());
    }
    return value;
}

class output_buffer {
public:
    void put(uint8_t c) {
        assert(used() < capacity());
        buffer_[used_++] = c;
    }

    void copy_match(int distance, int length) {
        assert(distance <= destination);
        assert(distance < 32768);
        assert(length >= 3 && length <= max_match_length);
        uint8_t* out = ptr() + used_;
        const uint8_t* in = out - distance;
        used_ += length;
        if (distance >= length) {
            memcpy(out, in, length);
        } else {
            do {
                *out++ = *in++;
            } while (--length);
        }
    }

    void add_used(int used) {
        assert(used <= avail());
        used_ += used;
    }

    int used() const {
        return used_;
    }

    int avail() const {
        return capacity() - used();
    }

    int capacity() const {
        return static_cast<int>(buffer_.size());
    }

    uint8_t* ptr() {
        return &buffer_[0];
    }

    void enlarge() {
        buffer_.resize(buffer_.size() < 32768 ? 32768 : buffer_.size() * 2);
    }

    auto steal_buffer() {
        buffer_.resize(used());
        return std::move(buffer_);
    }

private:
    std::vector<uint8_t> buffer_;
    int                  used_ = 0;
};

bool deflate_inner(output_buffer& output, bit_stream& bs, const huffman_tree& lit_len_tree, const huffman_tree& dist_tree)
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

    for (;;) {
        if (output.avail() < max_match_length) {
            output.enlarge();
        }

        // decode literal/length value from input stream
        const int value = decode(lit_len_tree, bs);
        if (value <= lit_max) {
            // if value < 256
            //    copy value (literal byte) to output stream
            output.put(static_cast<uint8_t>(value));
            //std::cout << "Lit " << litrep(value) << "\n";
        } else if (value == end_of_block) {
            // if value = end of block (256)
            //    break from loop
            return true;
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
            assert(len >= 3 && len <= max_match_length);

            int dist = decode(dist_tree, bs);
            const int dist_extra_bits = distance_extra_bits[dist];
            int dist_bytes = distance_length[dist];
            if (dist_extra_bits) {
                dist_bytes += bs.get_bits(dist_extra_bits);
            }

            if (dist_bytes > output.used()) {
                return false;
            }
            output.copy_match(dist_bytes, len);
        }
    }
}

bool read_dynamic_huffman_code_lengths(bit_stream& bs, std::vector<uint8_t>& cl2, int& number_of_lit_len_codes)
{
    // read representation of code trees
    const int hlit  = 257 + bs.get_bits(5);     // 5 Bits: HLIT, # of Literal/Length codes - 257 (257 - 286)
    const int hdist = 1 + bs.get_bits(5);       // 5 Bits: HDIST, # of Distance codes - 1        (1 - 32)
    const int hclen = 4 + bs.get_bits(4);       // 4 Bits: HCLEN, # of Code Length codes - 4     (4 - 19)

    constexpr int max_code_lengths = 19;
    uint8_t code_lengths[max_code_lengths]={0};
    for (int i = 0; i < hclen; ++i) {
        constexpr int alphabet_permute[max_code_lengths] = {
            16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
        };
        code_lengths[alphabet_permute[i]] = static_cast<uint8_t>(bs.get_bits(3));
    }

    huffman_tree cl_tree = make_huffman_tree(make_huffman_table(code_lengths, code_lengths+max_code_lengths), 7);

    number_of_lit_len_codes = hlit;
    cl2.resize(hlit + hdist);
    for (int i = 0; i < hlit + hdist;) {
        auto cl_val = static_cast<uint8_t>(decode(cl_tree, bs));
        int count = 1;
        if (cl_val <= 15) {
            // 0 - 15: Represent code lengths of 0 - 15
        } else if (cl_val == 16) {
            // 16: Copy the previous code length 3 - 6 times.
            if (i == 0) {
                return false;
            }
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
            if (cl_val != 18) {
                return false;
            }
            cl_val = 0;
            count = 11 + bs.get_bits(7);
        }
        if (cl_val > max_bits || static_cast<size_t>(count + i) > cl2.size()) {
            return false;
        }
        while (count--) {
            cl2[i++] = cl_val;
        }
    }

    return true;
}

std::vector<uint8_t> deflate(bit_stream& bs)
{
    auto invalid = [] () { assert(false); throw std::runtime_error("Invalid deflate stream"); };

    static const auto default_lit_len_tree = make_huffman_tree(make_default_huffman_table(), 9);
    static const auto default_dist_tree    = make_huffman_tree(make_default_huffman_len_table(), 5);

    output_buffer output;
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
        } else if (type == block_type::dynamic_huffman) {

            std::vector<uint8_t> cl2;
            int hlit;
            if (!read_dynamic_huffman_code_lengths(bs, cl2, hlit)) {
                invalid();
            }
            const auto lit_len_tree = make_huffman_tree(make_huffman_table(cl2.data(), cl2.data() + hlit), 9);
            const auto dist_tree    = make_huffman_tree(make_huffman_table(cl2.data() + hlit, cl2.data() + cl2.size()), 6);

            if (!deflate_inner(output, bs, lit_len_tree, dist_tree)) {
                invalid();
            }
        } else if (type == block_type::fixed_huffman) {
            if (!deflate_inner(output, bs, default_lit_len_tree, default_dist_tree)) {
                invalid();
            }
        } else {
            invalid();
        }
    }
    return output.steal_buffer();
}

} // namespace deflate
