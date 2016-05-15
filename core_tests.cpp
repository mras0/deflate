#include <string>
#include <iostream>
#include <algorithm>
#include <stdlib.h>

#include "crc.h"
#include "bit_stream.h"
#include "huffman_tree.h"
#include "huffman_table.h"
#include "deflate.h"

#define CHECK(expr) do { if (!(expr)) { std::cerr << #expr << std::endl; abort(); } } while (false)

using namespace deflate;

void test_crc32()
{
    const uint8_t d[14] = { 'L', 'i', 'n', 'e', ' ', '1', '\n', 'L', 'i', 'n', 'e', ' ', '2', '\n'};
    CHECK(0x87E4F545 == update_crc32(0, d, d+sizeof(d)));
}

void test_bit_stream()
{
    const uint8_t data[] = { 0x5a, 0xa5 }; // 01011010 10100101
    {
        bit_stream bs{data};
        CHECK(bs.potentially_available_bits() == 16);
        CHECK(bs.get_bits(16) == 0xa55a);
        CHECK(bs.available_bits() == 0);
    }
    {
        bit_stream bs{data};
        CHECK(bs.get_bits(8) == 0x5a);
        CHECK(bs.potentially_available_bits() == 8);
        CHECK(bs.get_bits(8) == 0xa5);
    }
    {
        bit_stream bs{data};
        CHECK(bs.get_bits(4) == 0xa);
        CHECK(bs.potentially_available_bits() == 12);
        CHECK(bs.get_bits(4) == 0x5);
    }
    {
        bit_stream bs{data};
        CHECK(bs.get_bits(2) == 0x2);
        CHECK(bs.get_bits(2) == 0x2);
        CHECK(bs.get_bits(2) == 0x1);
        CHECK(bs.get_bits(2) == 0x1);
    }
    {
        bit_stream bs{data};
        CHECK(bs.get_bit() == 0);
        CHECK(bs.get_bit() == 1);
        CHECK(bs.get_bit() == 0);
        CHECK(bs.available_bits() >= 5);
        CHECK(bs.potentially_available_bits() == 13);
        bs.ensure_bits(13);
        CHECK(bs.available_bits() >= 13);
        CHECK(bs.get_bit() == 1);
        CHECK(bs.get_bit() == 1);
        CHECK(bs.get_bit() == 0);
        CHECK(bs.get_bit() == 1);
        CHECK(bs.get_bit() == 0);
    }
}

std::ostream& operator<<(std::ostream& os, const huffman_code& c) {
    auto v = c.value;
    for (int i = c.len-1; i >= 0; --i) {
        os << ((v >> i) & 1);
    }
    return os;
}

void test_huffman_tree()
{
    auto te = [] (int len, int index) { return huffman_tree::table_entry{len, index}; };
    {
        constexpr huffman_code a_code{ 2, 0b00  };
        constexpr huffman_code b_code{ 1, 0b1   };
        constexpr huffman_code c_code{ 3, 0b011 };
        constexpr huffman_code d_code{ 3, 0b010 };
        huffman_tree t;
        t.add('A', a_code);
        t.add('B', b_code);
        t.add('C', c_code);
        t.add('D', d_code);
        CHECK(t.symbol(a_code) == 'A');
        CHECK(t.symbol(b_code) == 'B');
        CHECK(t.symbol(c_code) == 'C');
        CHECK(t.symbol(d_code) == 'D');
        CHECK(t.symbol_code('A') == a_code);
        CHECK(t.symbol_code('B') == b_code);
        CHECK(t.symbol_code('C') == c_code);
        CHECK(t.symbol_code('D') == d_code);
        t.make_tables(4);
        CHECK(t.next_from_bits(0b00, 4) == te(a_code.len, 'A'));
        CHECK(t.next_from_bits(0b1,  8) == te(b_code.len, 'B'));
        CHECK(t.next_from_bits(0b110, 4) == te(c_code.len, 'C'));
        CHECK(t.next_from_bits(0b010, 12) == te(d_code.len, 'D'));
    }
    {
        constexpr huffman_code a_code{ 2, 0b10  };
        constexpr huffman_code b_code{ 1, 0b0   };
        constexpr huffman_code c_code{ 3, 0b110 };
        constexpr huffman_code d_code{ 3, 0b111 };
        huffman_tree t;
        t.add('A', a_code);
        t.add('B', b_code);
        t.add('C', c_code);
        t.add('D', d_code);
        CHECK(t.symbol(a_code) == 'A');
        CHECK(t.symbol(b_code) == 'B');
        CHECK(t.symbol(c_code) == 'C');
        CHECK(t.symbol(d_code) == 'D');
        CHECK(t.symbol_code('A') == a_code);
        CHECK(t.symbol_code('B') == b_code);
        CHECK(t.symbol_code('C') == c_code);
        CHECK(t.symbol_code('D') == d_code);
        t.make_tables(2);
        CHECK(t.next_from_bits(0b01, 2) == te(2, 'A'));
        CHECK(t.next_from_bits(0b00, 2) == te(1, 'B'));
        CHECK(t.next_from_bits(0b11, 2).len == 2);
        CHECK(t.next_from_bits(0b11, 2).index >= huffman_tree::max_symbols);
        CHECK(t.branch(t.next_from_bits(0b11, 2).index-huffman_tree::max_symbols, false) == 'C');
        CHECK(t.branch(t.next_from_bits(0b11, 2).index-huffman_tree::max_symbols, true) == 'D');
    }
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    for (const auto& e : v) os << e << " ";
    return os;
}

void test_make_huffman_table()
{
    // Example from rfc1951 3.2.2
    std::vector<uint8_t> symbol_bit_lengths{3, 3, 3, 3, 3, 2, 4, 4};
    auto codes = make_huffman_table(symbol_bit_lengths);
    CHECK(codes == (std::vector<huffman_code>{
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
    for (int i = 0; i < huffman_tree::max_symbols; ++i) {
        auto expected = [i]() {
            if (i < 144) return huffman_code{8, static_cast<uint32_t>(0b00110000  + (i - 0))};
            if (i < 256) return huffman_code{9, static_cast<uint32_t>(0b110010000 + (i - 144))};
            if (i < 280) return huffman_code{7, static_cast<uint32_t>(0b0000000   + (i - 256))};
            assert(i < huffman_tree::max_symbols);
            return huffman_code{8, static_cast<uint32_t>(0b11000000  + (i - 280))};
        }();
        CHECK(cs[i] == expected);
    }
    huffman_tree t2;
    for (int i = 0; i < huffman_tree::max_symbols; ++i) {
        t2.add(i, cs[i]);
    }
}


void test_deflate()
{
    const uint8_t deflate_input1[13] = {0xf3, 0xc9, 0xcc, 0x4b, 0x55, 0x30, 0xe4, 0xf2, 0x01, 0x51, 0x46, 0x5c, 0x00};
    const uint8_t deflate_input2[12] = {0xf3, 0xc9, 0xcc, 0x4b, 0x55, 0x30, 0xe4, 0x02, 0x53, 0x46, 0x5c, 0x00};
    const std::vector<uint8_t> expected_output{ 'L', 'i', 'n', 'e', ' ', '1', '\n', 'L', 'i', 'n', 'e', ' ', '2', '\n'};
    bit_stream bs1{deflate_input1};
    CHECK(deflate::deflate(bs1) == expected_output);
    bit_stream bs2{deflate_input2};
    CHECK(deflate::deflate(bs2) == expected_output);
}

int main()
{
    try {
        test_crc32();
        test_bit_stream();
        test_huffman_tree();
        test_make_huffman_table();
        test_deflate();
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
