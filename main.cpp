#include <iostream>
#include <iomanip>
#include <string>
#include <stdint.h>
#include <cassert>
#include <array>

constexpr uint32_t crc32_poly = 0xedb88320; // 0x04C11DB7 reversed

constexpr uint32_t crc32_one_bit(uint32_t c)
{
    return c & 1 ? crc32_poly ^ (c >> 1) : c >> 1;
}

constexpr uint32_t crc32_one_byte(uint32_t c)
{
    return crc32_one_bit(crc32_one_bit(crc32_one_bit(crc32_one_bit(crc32_one_bit(crc32_one_bit(crc32_one_bit(crc32_one_bit(c))))))));
}

#define USE_CRC32_TABLE
#ifdef USE_CRC32_TABLE
template<size_t... is>
constexpr std::array<uint32_t, 256> make_crc32_table(std::index_sequence<is...>)
{
    return { crc32_one_byte(is)... };
}

constexpr std::array<uint32_t, 256> crc32_table = make_crc32_table(std::make_index_sequence<256>());
#endif

uint32_t crc32(uint32_t crc, const uint8_t* beg, const uint8_t* end)
{
    crc ^= ~0U;
    for (auto p = beg; p != end; ++p) {
#ifdef USE_CRC32_TABLE
        crc = crc32_table[(crc ^ *p) & 0xff] ^ (crc >> 8);
#else
        crc = crc32_one_byte(crc ^ *p);
#endif
    }
    crc ^= ~0U;
    return crc;
}

void test_crc32()
{
    const uint8_t d[14] = { 'L', 'i', 'n', 'e', ' ', '1', '\n', 'L', 'i', 'n', 'e', ' ', '2', '\n'};
    assert(0x87E4F545 == crc32(0, d, d+sizeof(d)));
}

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

std::string litrep(int value)
{
    std::string res;
    value = static_cast<unsigned char>(value);
    if (value >= 32 && value < 127) {
        res += static_cast<char>(value);
    }
    else {
        const char* hex = "0123456789abcdef";
        res += '\\';
        res += 'x';
        res += hex[value>>4];
        res += hex[value&15];
    }
    return res;
}

static constexpr int num_symbols = 288;
static constexpr int max_bits    = 15;

class huffman_tree {
public:
    explicit huffman_tree() {
    }

    int& branch(int index, bool right) {
        assert(index < num_nodes);
        auto& n = nodes[index];
        return right ? n.right : n.left;
    }

    const int& branch(int index, bool right) const {
        return const_cast<huffman_tree&>(*this).branch(index, right);
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

    void show(std::ostream& os) const {
        for (int i = 0; i < num_nodes; ++i) {
            os << std::setw(3) << i << " ";
            if (nodes[i].left >= num_symbols) os << "*" << (nodes[i].left - num_symbols);
            else os << litrep(nodes[i].left);
            os << " ";
            if (nodes[i].right >= num_symbols) os << "*" << (nodes[i].right - num_symbols);
            else os << litrep(nodes[i].right);
            os << "\n";
        }
    }

private:
    static constexpr int max_nodes          = num_symbols;
    static constexpr int invalid_edge_value = num_symbols + max_nodes;

    struct node {
        int left;
        int right;
    };

    int num_nodes = 0;
    node nodes[max_nodes];

    int alloc_node() {
        assert(num_nodes < max_nodes);
        nodes[num_nodes].left  = invalid_edge_value;
        nodes[num_nodes].right = invalid_edge_value;
        return (num_nodes++) + num_symbols;
    }

    static code bit_added(const code& c, uint8_t bit) {
        assert(c.len < UINT8_MAX - 1);
        assert(bit == 0 || bit == 1);
        return { static_cast<uint8_t>(c.len + 1), (c.value << 1) | static_cast<uint32_t>(bit) };
    }

    static bool consume_bit(code& c) {
        assert(c.valid());
        const auto mask = 1 << (c.len - 1);
        const bool ret  = (c.value & mask) != 0;
        c.value &= ~mask;
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
        constexpr code c_code{ 3, 0b011 };
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
        constexpr code a_code{ 2, 0b10  };
        constexpr code b_code{ 1, 0b0   };
        constexpr code c_code{ 3, 0b110 };
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

std::vector<code> make_huffman_table(const uint8_t* symbol_bit_lengths, const uint8_t* symbol_bit_lengths_end)
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

    std::vector<::code> codes(symbol_bit_lengths_end - symbol_bit_lengths);
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

std::vector<code> make_huffman_table(const std::vector<uint8_t>& symbol_bit_lengths)
{
    return make_huffman_table(symbol_bit_lengths.data(), symbol_bit_lengths.data() + symbol_bit_lengths.size());
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

auto make_default_huffman_len_table()
{
    // Distance codes 0-31 are represented by (fixed-length) 5-bit codes
    std::vector<uint8_t> length_bit_lengths(32, 5);
    return make_huffman_table(length_bit_lengths);
}

huffman_tree make_huffman_tree(const std::vector<code>& codes)
{
    assert(codes.size() <= num_symbols);
    huffman_tree t;
    for (int i = 0; i < (int)codes.size(); ++i) {
        if (codes[i].len > 0) {
            t.add(i, codes[i]);
        }
    }
    return t;
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

    huffman_tree t;
    for (int i = 0; i < (int)symbol_bit_lengths.size(); ++i) {
        t.add('A' + i, codes[i]);
    }

    auto cs = make_default_huffman_table();
    for (int i = 0; i < num_symbols; ++i) {
        auto expected = [i]() {
            if (i < 144) return code{8, static_cast<uint32_t>(0b00110000  + (i - 0))};
            if (i < 256) return code{9, static_cast<uint32_t>(0b110010000 + (i - 144))};
            if (i < 280) return code{7, static_cast<uint32_t>(0b0000000   + (i - 256))};
            assert(i < num_symbols);
            return code{8, static_cast<uint32_t>(0b11000000  + (i - 280))};
        }();
        assert(cs[i] == expected);
    }
    huffman_tree t2;
    for (int i = 0; i < num_symbols; ++i) {
        t2.add(i, cs[i]);
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

int decode(const huffman_tree& t, bit_stream& bs)
{
    int value = 0;
    for (;;) {
        value = t.branch(value, !!bs.get_bit());
        if (value < num_symbols) {
            break;
        }
        value -= num_symbols;
    }
    return value;
}

std::vector<uint8_t> deflate(bit_stream& bs)
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

                huffman_tree cl_tree = make_huffman_tree(make_huffman_table(code_lengths, code_lengths+max_code_lengths));

                std::vector<uint8_t> cl2(hlit + hdist);
                for (int i = 0; i < hlit + hdist;) {
                    auto cl_val = static_cast<uint8_t>(decode(cl_tree, bs));
                    int count = 1;
                    if (cl_val <= max_bits) {
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
                    if (count + i > cl2.size()) invalid();
                    while (count--) {
                        //std::cout << (i > hlit ? "dist" : "lit") <<  "[" << (i > hlit ? i - hlit : i) << "] = " << (int)cl_val << std::endl;
                        cl2[i++] = cl_val;
                    }
                }

                lit_len_tree = make_huffman_tree(make_huffman_table(cl2.data(), cl2.data() + hlit));
                dist_tree    = make_huffman_tree(make_huffman_table(cl2.data() + hlit, cl2.data() + hlit + hdist));
            } else {
                assert(type == block_type::fixed_huffman);
                lit_len_tree = make_huffman_tree(make_default_huffman_table());
                dist_tree    = make_huffman_tree(make_default_huffman_len_table());
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
                    if (dist_bytes > output.size()) invalid();
                    auto src = output.size() - dist_bytes;
                    for (int i = 0; i < len; ++i) {
                        //std::cout << litrep(output[src]);
                        output.push_back(output[src++]);
                    }
                    //std::cout << "\n";
                }
            }
        }
    }
    return output;
}

void test_deflate()
{
    const uint8_t deflate_input1[13] = {0xf3, 0xc9, 0xcc, 0x4b, 0x55, 0x30, 0xe4, 0xf2, 0x01, 0x51, 0x46, 0x5c, 0x00};
    const uint8_t deflate_input2[12] = {0xf3, 0xc9, 0xcc, 0x4b, 0x55, 0x30, 0xe4, 0x02, 0x53, 0x46, 0x5c, 0x00};
    const std::vector<uint8_t> expected_output{ 'L', 'i', 'n', 'e', ' ', '1', '\n', 'L', 'i', 'n', 'e', ' ', '2', '\n'};
    bit_stream bs1{deflate_input1};
    assert(deflate(bs1) == expected_output);
    bit_stream bs2{deflate_input2};
    assert(deflate(bs2) == expected_output);
}

#include <fstream>
void gunzip(const std::string& filename)
{
    std::ifstream in(filename, std::ios_base::binary);
    if (!in) throw std::runtime_error(filename + " not found");
    in.seekg(0, std::ios_base::end);
    const auto file_size = static_cast<int>(in.tellg());
    if (file_size < 18) throw std::runtime_error(filename + " is too small to be a gzip file");
    in.seekg(0, std::ios_base::beg);
    std::vector<uint8_t> input(file_size);
    in.read(reinterpret_cast<char*>(&input[0]), input.size());


    auto invalid = [&filename] () { throw std::runtime_error(filename + " is not a valid gzip file"); };

    // +---+---+---+---+---+---+---+---+---+---+
    // |ID1|ID2|CM |FLG|     MTIME      |XFL|OS|
    // +---+---+---+---+---+---+---+---+---+---+
    if (input[0] != 0x1f || // ID1
        input[1] != 0x8b || // ID2
        input[2] != 8) {    // CM
        invalid();
    }
    enum { FTEXT=1, FHCRC=2, FEXTRA=4, FNAME=8, FCOMMENT=16 };
    const auto flg = input[3];
    int pos = 10;
    if (flg & FEXTRA) {
        const int xlen = input[pos] | (input[pos+1] << 8);
        // TODO: Check whether it's valid to skip this many bytes
        pos += 2 + xlen;
    }
    auto get_zero_terminated_string = [&] () -> std::string {
        std::string res;
        while (pos < file_size) {
            const char c = input[pos++];
            if (!c) return res;
            res += c;
        }
        invalid();
        return {};
    };
    if (flg & FNAME) {
        std::cout << "Filename: " << get_zero_terminated_string() << std::endl;
    }
    if (flg & FCOMMENT) {
        std::cout << "Comment: " << get_zero_terminated_string() << std::endl;
    }
    if (flg & FHCRC) {
        pos += 2; // skip 16-bit header CRC
    }
    assert(file_size >= pos + 8);

    int f = file_size - 8;
    const auto crc32 = static_cast<uint32_t>(input[f + 0]) |
        (static_cast<uint32_t>(input[f + 1]) << 8) |
        (static_cast<uint32_t>(input[f + 2]) << 16) |
        (static_cast<uint32_t>(input[f + 3]) << 24);

    const auto isize = static_cast<uint32_t>(input[f + 4]) |
        (static_cast<uint32_t>(input[f + 5]) << 8) |
        (static_cast<uint32_t>(input[f + 6]) << 16) |
        (static_cast<uint32_t>(input[f + 7]) << 24);

    bit_stream bs(input.data() + pos, input.data() + file_size - 8);
    const auto output = deflate(bs);
    if (output.size() != isize) {
        invalid();
    }
    std::cout.write(reinterpret_cast<const char*>(output.data()), output.size());

    if (::crc32(0, output.data(), output.data() + output.size()) != crc32) {
        invalid();
    }
}

int main()
{
    try {
        test_crc32();
        test_bit_stream();
        test_huffman_tree();
        test_make_huffman_table();
        test_deflate();
        gunzip("../CMakeLists.txt.gz");
        gunzip("../main.cpp.gz");
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
