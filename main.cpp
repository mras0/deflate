#include <iostream>
#include <iomanip>
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

    bool valid() const {
        return len > 0 && len <= 8 * sizeof(value) && (value >> len) == 0;

    }
};
bool operator==(const code& l, const code& r) { return l.len == r.len && l.value == r.value; }
bool operator!=(const code& l, const code& r) { return !(l == r); }
std::ostream& operator<<(std::ostream& os, const code& c) {
    auto v = c.value;
    for (int i = c.len-1; i >= 0; --i) {
        os << ((v >> i) & 1);
    }
    return os;
}

static constexpr int num_symbols = 288;

class huffman_tree {
public:
    explicit huffman_tree() {
    }

    void add(int symbol, const code& symbol_code) {
        assert(symbol_code.valid());
        auto c = symbol_code;

        if (num_nodes == 0) {
            alloc_node();
        }

        int index = 0;
        // Create non-leaf nodes
        while (c.len > 1) {
            auto& edge = branch(index, consume_bit(c));
            if (edge == invalid_edge_value) {
                edge = alloc_node();
            }
            assert(edge >= num_symbols);
            index = edge - num_symbols;
        }
        auto& edge = branch(index, consume_bit(c));
        assert(edge == invalid_edge_value);
        edge = symbol;
    }

    int symbol(const code& symbol_code) const {
        assert(symbol_code.valid());
        auto c = symbol_code;
        int index = 0;
        while (c.len > 1) {
            index = branch(index, consume_bit(c));
            assert(index >= num_symbols);
            index -= num_symbols;
        }

        const auto ret = branch(index, consume_bit(c));
        assert(ret < num_symbols);
        return ret;
    }

    code symbol_code(int symbol) const {
        code c{};
        coder(symbol, 0, c);
        assert(c.valid());
        return c;
    }

private:
    static constexpr int max_nodes          = 32;
    static constexpr int invalid_edge_value = num_symbols + max_nodes;

    struct node {
        int left;
        int right;
    };

    int num_nodes = 0;
    node nodes[max_nodes];

    int& branch(int index, bool right) {
        assert(index < num_nodes);
        auto& n = nodes[index];
        return right ? n.right : n.left;
    }
    const int& branch(int index, bool right) const {
        return const_cast<huffman_tree&>(*this).branch(index, right);
    }

    int alloc_node() {
        assert(num_nodes < max_nodes);
        nodes[num_nodes].left  = invalid_edge_value;
        nodes[num_nodes].right = invalid_edge_value;
        return (num_nodes++) + num_symbols;
    }

    static code bit_added(const code& c, uint8_t bit) {
        assert(c.len < UINT8_MAX - 1);
        assert(bit == 0 || bit == 1);
        return { static_cast<uint8_t>(c.len + 1), c.value | (static_cast<uint32_t>(bit) << c.len) };
    }

    static bool consume_bit(code& c) {
        assert(c.valid());
        const bool ret = (c.value & 1) != 0;
        c.value >>= 1;
        c.len--;
        return ret;
    }

    bool coder(int symbol, int index, code& c) const {
        const auto orig = c;

        assert(index < num_nodes);
        auto& n = nodes[index];
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
};

void test_huffman_tree()
{
    {
        constexpr code a_code{ 2, 0b00  };
        constexpr code b_code{ 1, 0b1   };
        constexpr code c_code{ 3, 0b110 };
        constexpr code d_code{ 3, 0b010 };
        huffman_tree t;
        t.add('A', a_code);
        t.add('B', b_code);
        t.add('C', c_code);
        t.add('D', d_code);
        assert(t.symbol(a_code) == 'A');
        assert(t.symbol(b_code) == 'B');
        assert(t.symbol(c_code) == 'C');
        assert(t.symbol(d_code) == 'D');
        assert(t.symbol_code('A') == a_code);
        assert(t.symbol_code('B') == b_code);
        assert(t.symbol_code('C') == c_code);
        assert(t.symbol_code('D') == d_code);
    }
    {
        constexpr code a_code{ 2, 0b01  };
        constexpr code b_code{ 1, 0b0   };
        constexpr code c_code{ 3, 0b011 };
        constexpr code d_code{ 3, 0b111 };
        huffman_tree t;
        t.add('A', a_code);
        t.add('B', b_code);
        t.add('C', c_code);
        t.add('D', d_code);
        assert(t.symbol(a_code) == 'A');
        assert(t.symbol(b_code) == 'B');
        assert(t.symbol(c_code) == 'C');
        assert(t.symbol(d_code) == 'D');
        assert(t.symbol_code('A') == a_code);
        assert(t.symbol_code('B') == b_code);
        assert(t.symbol_code('C') == c_code);
        assert(t.symbol_code('D') == d_code);
    }
}

#include <vector>
#include <algorithm>

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    for (const auto& e : v) os << e << " ";
    return os;
}

std::vector<code> make_huffman_table(const std::vector<uint8_t>& symbol_bit_lengths)
{
    const auto max_bit_length = *std::max_element(begin(symbol_bit_lengths), end(symbol_bit_lengths));
    // Count the number of codes for each code length. Let bl_count[N] be the number of codes of length N, N >= 1.
    std::vector<int> bl_count(max_bit_length+1);
    for (const auto& bl : symbol_bit_lengths) {
        ++bl_count[bl];
    }

    //  Find the numerical value of the smallest code for each code length
    int code = 0;
    assert(bl_count[0] == 0);
    std::vector<uint32_t> next_code(max_bit_length+1);
    for (int bits = 1; bits <= max_bit_length; bits++) {
        code = (code + bl_count[bits-1]) << 1;
        next_code[bits] = code;
    }

    // Assign numerical values to all codes, using consecutive values for all codes of the same length with the base
    // values determined at step 2. Codes that are never used (which have a bit length of zero) must not be assigned
    // a value.

    std::vector<::code> codes(symbol_bit_lengths.size());
    for (int n = 0; n < (int)codes.size(); n++) {
        const auto len = symbol_bit_lengths[n];//len = tree[n].Len;
        if (len != 0) {
            //tree[n].Code = next_code[len];
            codes[n] = ::code{len, next_code[len]};
            next_code[len]++;
        }
    }
    return codes;
}

auto make_default_huffman_table()
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
    std::vector<uint8_t> sbl(num_symbols);
    int i = 0;
    while (i < 144) sbl[i++] = 8;
    while (i < 256) sbl[i++] = 9;
    while (i < 280) sbl[i++] = 7;
    while (i < num_symbols) sbl[i++] = 8;
    return make_huffman_table(sbl);
}

void test_make_huffman_table()
{
    // Example from rfc1951 3.2.2
    std::vector<uint8_t> symbol_bit_lengths{3, 3, 3, 3, 3, 2, 4, 4};
    auto codes = make_huffman_table(symbol_bit_lengths);
    assert(codes == (std::vector<code>{
        { 3,  0b010 },
        { 3,  0b011 },
        { 3,  0b100 },
        { 3,  0b101 },
        { 3,  0b110 },
        { 2,   0b00 },
        { 4, 0b1110 },
        { 4, 0b1111 }}));

   auto cs = make_default_huffman_table();
   for (int i = 0; i < num_symbols; ++i) {
       auto expected = [i] () {
           if (i < 144) return code{8, static_cast<uint32_t>(0b00110000  + (i - 0))};
           if (i < 256) return code{9, static_cast<uint32_t>(0b110010000 + (i - 144))};
           if (i < 280) return code{7, static_cast<uint32_t>(0b0000000   + (i - 256))};
           assert(i < num_symbols);
           return code{8, static_cast<uint32_t>(0b11000000  + (i - 280))};
       }();
       const auto& c = cs[i];
       assert(c == expected);
   }
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
    test_make_huffman_table();
    const uint8_t deflate_input1[13] = {0xf3, 0xc9, 0xcc, 0x4b, 0x55, 0x30, 0xe4, 0xf2, 0x01, 0x51, 0x46, 0x5c, 0x00};
    const uint8_t deflate_input2[12] = {0xf3, 0xc9, 0xcc, 0x4b, 0x55, 0x30, 0xe4, 0x02, 0x53, 0x46, 0x5c, 0x00};
    const uint8_t expected_output[14] = { 'L', 'i', 'n', 'e', ' ', '1', '\n', 'L', 'i', 'n', 'e', ' ', '2', '\n'};
    bit_stream bs{deflate_input1};
    //deflate(bs);
}
