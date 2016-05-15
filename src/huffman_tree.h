#ifndef DEFLATE_HUFFMAN_TREE_H
#define DEFLATE_HUFFMAN_TREE_H

#include <stdint.h>
#include <cassert>
#include <iosfwd>
#include <vector>

#include "huffman_code.h"

namespace deflate {

class huffman_tree {
public:
    static constexpr int max_symbols = 288;

    explicit huffman_tree() {
    }

    // branch return value meaning:
    //  0..max_symbols-1 symbol
    //  max_symbols..    internal node (next index = return value - max_symbols)
    uint16_t& branch(int index, bool right) {
        assert(index < num_nodes);
        auto& n = nodes[index];
        return right ? n.right : n.left;
    }

    const uint16_t& branch(int index, bool right) const {
        return const_cast<huffman_tree&>(*this).branch(index, right);
    }

    void add(int symbol, const huffman_code& symbol_code) {
        assert(symbol >= 0 && symbol < max_symbols);
        assert(symbol_code.valid());
        assert(symbol_code.len <= max_bits);
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
            assert(edge >= max_symbols);
            index = edge - max_symbols;
        }
        auto& edge = branch(index, consume_bit(c));
        assert(edge == invalid_edge_value);
        edge = static_cast<uint16_t>(symbol);
    }

    int symbol(const huffman_code& symbol_code) const {
        assert(symbol_code.valid());
        auto c = symbol_code;
        int index = 0;
        while (c.len > 1) {
            index = branch(index, consume_bit(c));
            assert(index >= max_symbols);
            index -= max_symbols;
        }

        const auto ret = branch(index, consume_bit(c));
        assert(ret < max_symbols);
        return ret;
    }

    huffman_code symbol_code(int symbol) const {
        huffman_code c{};
        coder(symbol, 0, c);
        assert(c.valid());
        return c;
    }

    void output_graph(std::ostream& os) const;

    class table_entry {
    public:
        explicit table_entry() : repr_(0) {
        }        
        explicit table_entry(int len, int index) : repr_(static_cast<uint16_t>((len<<12)|index)) {
            assert(len > 0 && len <= max_bits);
            assert(index < invalid_edge_value);
        }
        explicit table_entry(uint16_t repr) : repr_(repr) {
        }
        uint8_t  len() const { return static_cast<uint8_t>(repr_ >> 12); }
        uint16_t index() const { return static_cast<uint16_t>(repr_ & 0xfff); }

        bool operator==(const table_entry& rhs) const { return repr_ == rhs.repr_; }

    private:
        uint16_t repr_;
    };

    table_entry next_from_bits(uint32_t bits, int num_bits) const {
        assert(table_bits());
        assert(num_bits >= table_bits()); (void)num_bits;
        return table_entry { table[bits & ((1 << table_bits_) -1)] };
    }

    int table_bits() const {
        return table_bits_;
    }

    void make_tables(int num_table_bits) {
        assert(num_table_bits > 0 && num_table_bits <= max_table_bits);
        assert(num_nodes > 0);
        const int table_size = 1 << num_table_bits;
        table_bits_ = num_table_bits;
        for (int i = 0; i < table_size; ++i) {
            int index = max_symbols;
            int len = 0;
            int val = i;
            while (len < table_bits_ && index >= max_symbols) {
                ++len;
                index = branch(index - max_symbols, val & 1);
                val >>= 1;
            }
            table[i] = table_entry{len, index};
        }
    }

private:
    static constexpr int max_nodes          = max_symbols;
    static constexpr int invalid_edge_value = max_symbols + max_nodes;
    static constexpr int max_table_bits     = 9;

    struct node {
        uint16_t left;
        uint16_t right;
    };

    int num_nodes = 0;
    node nodes[max_nodes];

    int table_bits_ = 0;
    table_entry table[1<<max_table_bits];

    uint16_t alloc_node() {
        assert(num_nodes < max_nodes);
        nodes[num_nodes].left  = invalid_edge_value;
        nodes[num_nodes].right = invalid_edge_value;
        return static_cast<uint16_t>((num_nodes++) + max_symbols);
    }

    static huffman_code bit_added(const huffman_code& c, uint8_t bit) {
        assert(c.len <= max_bits);
        assert(bit == 0 || bit == 1);
        return { static_cast<uint8_t>(c.len + 1), static_cast<uint16_t>((c.value << 1) | bit) };
    }

    static bool consume_bit(huffman_code& c) {
        assert(c.valid());
        const auto mask = 1 << (c.len - 1);
        const bool ret  = (c.value & mask) != 0;
        c.value &= ~mask;
        c.len--;
        return ret;
    }

    bool coder(int symbol, int index, huffman_code& c) const {
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
        if (n.left >= max_symbols && coder(symbol, n.left - max_symbols, c)) {
            return true;
        }
        c = r;
        if (n.right >= max_symbols && coder(symbol, n.right - max_symbols, c)) {
            return true;
        }

        c = orig;
        return false;
    }
};

huffman_tree make_huffman_tree(const std::vector<huffman_code>& codes, int table_bits);

} // namespace deflate

#endif
