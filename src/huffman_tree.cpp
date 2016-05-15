#include "huffman_tree.h"
#include <ostream>
#include <string>

namespace deflate {

void huffman_tree::output_graph(std::ostream& os) const
{
    os << "digraph G {\n";

    auto wr_val = [&](int val, const char* lab) {
        std::string ltext = " [label=\"" + std::string(lab) + "\"]";
        if (val >= max_symbols) os << "node" << (val - max_symbols) << ltext;
        else if ((val >= 'A' && val <= 'Z') || (val >= 'a' && val <= 'z')) os << (char)val << ltext;
        else os << "val" << val <<  ltext << "\nval" << val << "[label=\"" << val << "\"]";
    };

    for (int i = 0; i < num_nodes; ++i) {
        os << "node" << i << " [label=\"\"]\n";

        os << "node" << i << " -> ";
        wr_val(nodes[i].left, "0");
        os << "\n";

        os << "node" << i << " -> ";
        wr_val(nodes[i].right, "1");
        os << "\n";
    }
    os << "}\n";
}

huffman_tree make_huffman_tree(const std::vector<huffman_code>& codes, int table_bits)
{
    assert(codes.size() <= huffman_tree::max_symbols);
    huffman_tree t;
    for (int i = 0; i < (int)codes.size(); ++i) {
        if (codes[i].len > 0) {
            t.add(i, codes[i]);
        }
    }
    t.make_tables(table_bits);
    return t;
}

} // namespace deflate