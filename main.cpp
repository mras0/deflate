#include <iostream>
#include <stdint.h>
#include <cassert>

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
    template<int size>
    explicit bit_stream(const uint8_t (&data)[size]) : data_(data), len_(size) {
    }

    uint8_t get_bit() {
        return static_cast<uint8_t>(get_bits(1));
    }

    uint32_t get_bits(int num_bits) {
        assert(num_bits > 0 && num_bits < 24);
        while (avail_ < num_bits) {
            bits_   = (static_cast<uint32_t>(get_byte()) << avail_) | bits_;
            avail_ += 8;
        }
        const auto res = bits_ & ((1 << num_bits) - 1);
        avail_ -= num_bits;
        bits_ >>= num_bits;
        return res;
    }

private:
    const uint8_t*  data_;
    int             len_;
    int             pos_   = 0;
    uint32_t        bits_  = 0;
    int             avail_ = 0;

    uint8_t get_byte() {
        if (pos_ < len_) {
            return data_[pos_++];
        }
        assert(false);
        return 0;
    }
};

void test_bit_stream()
{
    const uint8_t data[] = { 0x5a, 0xa5 }; // 01011010 10100101
    {
        bit_stream bs{data};
        assert(bs.get_bits(16) == 0xa55a);
    }
    {
        bit_stream bs{data};
        assert(bs.get_bits(8) == 0x5a);
        assert(bs.get_bits(8) == 0xa5);
    }
    {
        bit_stream bs{data};
        assert(bs.get_bits(4) == 0xa);
        assert(bs.get_bits(4) == 0x5);
    }
    {
        bit_stream bs{data};
        assert(bs.get_bits(2) == 0x2);
        assert(bs.get_bits(2) == 0x2);
        assert(bs.get_bits(2) == 0x1);
        assert(bs.get_bits(2) == 0x1);
    }
    {
        bit_stream bs{data};
        assert(bs.get_bit() == 0);
        assert(bs.get_bit() == 1);
        assert(bs.get_bit() == 0);
        assert(bs.get_bit() == 1);
        assert(bs.get_bit() == 1);
        assert(bs.get_bit() == 0);
        assert(bs.get_bit() == 1);
        assert(bs.get_bit() == 0);
    }
}

struct code {
    uint8_t  len;
    uint32_t value;
};
bool operator==(const code& l, const code& r) { return l.len == r.len && l.value == r.value; }
bool operator!=(const code& l, const code& r) { return !(l == r); }
std::ostream& operator<<(std::ostream& os, const code& c) {
    auto v = c.value;
    for (int i = 0; i < c.len; ++i) {
        os << (v & 1);
        v >>= 1;
    }
    return os;
}

code bit_added(const code& c, uint8_t bit) {
    assert(c.len < UINT8_MAX - 1);
    assert(bit == 0 || bit == 1);
    return { static_cast<uint8_t>(c.len + 1), (c.value << 1) | bit };
}

class huffman_tree {
public:
    explicit huffman_tree() {
        nodes_[0].left   = num_symbols+1;
        nodes_[0].right  = 'B';
        nodes_[1].left   = 'A';
        nodes_[1].right  = num_symbols+2;
        nodes_[2].left   = 'D';
        nodes_[2].right  = 'C';
    }

    int symbol(const code& c) const {
        assert(c.len > 0 && c.len <= 8 * sizeof(c.value));
        assert((c.value >> c.len) == 0);
        int index = 0;
        for (int bit = c.len - 1; bit > 0; --bit) {
            index = node_value(index, ((c.value >> bit) & 1) != 0);
            assert(index >= num_symbols);
            index -= num_symbols;
        }
        const auto ret = node_value(index, (c.value & 1) != 0);
        assert(ret < num_symbols);
        return ret;
    }

    code symbol_code(int symbol) const {
        code c{};
        coder(symbol, 0, c);
        return c;
    }

private:
    static constexpr int num_symbols = 256;
    static constexpr int max_nodes   = 3;//285-256;

    struct node {
        int left;
        int right;
    };

    node nodes_[max_nodes];

    bool coder(int symbol, int index, code& c) const {
        const auto orig = c;

        assert(index < max_nodes);
        auto& n = nodes_[index];
        const auto l = bit_added(c, 0);
        const auto r = bit_added(c, 1);

        // Simple matches
        if (n.left == symbol) {
            c = l;
            return true;
        } else if (n.right == symbol) {
            c = r;
            return true;
        }

        c = l;
        if (n.left >= num_symbols && coder(symbol, n.left - num_symbols, c)) {
            return true;
        }
        c = r;
        if (n.right >= num_symbols && coder(symbol, n.right - num_symbols, c)) {
            return true;
        }

        c = orig;
        return false;
    }

    int node_value(int index, bool right) const {
        assert(index < max_nodes);
        auto& n = nodes_[index];
        return right ? n.right : n.left;
    }
};

void test_huffman_tree()
{
    constexpr code a_code{ 2, 0b00  };
    constexpr code b_code{ 1, 0b1   };
    constexpr code c_code{ 3, 0b011 };
    constexpr code d_code{ 3, 0b010 };
    huffman_tree t;
    assert(t.symbol(a_code) == 'A');
    assert(t.symbol(b_code) == 'B');
    assert(t.symbol(c_code) == 'C');
    assert(t.symbol(d_code) == 'D');
    assert(t.symbol_code('A') == a_code);
    assert(t.symbol_code('B') == b_code);
    assert(t.symbol_code('C') == c_code);
    assert(t.symbol_code('D') == d_code);
#if 0
        A     
        B     0b1
        C     0b011
        D     0b010
#endif
}

enum class block_type { uncompressed, fixed_huffman, dynamic_huffman, reserved };
std::ostream& operator<<(std::ostream& os, block_type t) {
    switch (t) {
    case block_type::uncompressed:      return os << "uncompressed";
    case block_type::fixed_huffman:     return os << "fixed_huffman";
    case block_type::dynamic_huffman:   return os << "dynamic_huffman";
    case block_type::reserved:          return os << "reserved";
    }
    assert(false);
    return os << "block_type{" << static_cast<int>(t) << "}";
}

void deflate(bit_stream& bs)
{
    constexpr int min_match_distance_bytes          = 1;
    constexpr int max_match_distance_bytes          = 32768;
    constexpr int min_match_length_bytes            = 3;
    constexpr int max_match_length_bytes            = 258;
    enum alphabet {
        lit_min      = 0,
        lit_max      = 255,
        end_of_block = 256,
        len_min      = 257,
        len_max      = 285,
    };

    bool last_block = false;
    while (!last_block) {
        // Block header bits
        last_block = !!bs.get_bit();
        const auto type  = static_cast<block_type>(bs.get_bits(2));
        std::cout << "final " << last_block << " type " << type << std::endl;
        assert(type != block_type::reserved);

        if (type == block_type::uncompressed) {
            assert(false);
            // skip any remaining bits in current partially
            //     processed byte
            //     read LEN and NLEN (each 16-bits)
            //     copy LEN bytes of data to output
        } else {
            if (type == block_type::dynamic_huffman) {
                assert(false);
                // read representation of code trees
            }

            bool end_of_block = false;
            while (!end_of_block) {
                // decode literal/length value from input stream
                for (int i = 0; i < 9; ++i) std::cout << (int)bs.get_bit();
                std::cout << std::endl;
                
                assert(false);
                // if value < 256
                //    copy value (literal byte) to output stream
                // otherwise
                //    if value = end of block (256)
                //       break from loop
                //    otherwise (value = 257..285)
                //       decode distance from input stream
                // 
                //       move backwards distance bytes in the output
                //       stream, and copy length bytes from this
                //       position to the output stream.
            }
        }
    }
}

int main()
{
    test_bit_stream();
    test_huffman_tree();
    const uint8_t deflate_input1[13] = {0xf3, 0xc9, 0xcc, 0x4b, 0x55, 0x30, 0xe4, 0xf2, 0x01, 0x51, 0x46, 0x5c, 0x00};
    const uint8_t deflate_input2[12] = {0xf3, 0xc9, 0xcc, 0x4b, 0x55, 0x30, 0xe4, 0x02, 0x53, 0x46, 0x5c, 0x00};
    const uint8_t expected_output[14] = { 'L', 'i', 'n', 'e', ' ', '1', '\n', 'L', 'i', 'n', 'e', ' ', '2', '\n'};
    bit_stream bs{deflate_input1};
    deflate(bs);
}
