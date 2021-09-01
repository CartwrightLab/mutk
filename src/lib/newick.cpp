/*
# Copyright (c) 2020 Reed A. Cartwright <reed@cartwright.ht>
#
# This file is part of the Ultimate Source Code Project.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
*/

#include <doctest/doctest.h>

#include <mutk/detail/graph.hpp>
#include <mutk/utility.hpp>

#include <memory>

#include <boost/spirit/home/x3.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/sequence.hpp>

using mutk::utility::percent_decode;

namespace newick {
// http://evolution.genetics.washington.edu/phylip/newick_doc.html
namespace x3 = boost::spirit::x3;
namespace ascii = boost::spirit::x3::ascii;

struct node_data_t {
    node_data_t() = default;
    node_data_t(std::string l, float f) :
        label{l}, length{f} , max_length{f}, total_length{f}
    {}

    std::string label; // node label
    float length; // length from parent node
    float max_length;
    float total_length;
    std::size_t parent{0}; // index of parent
};

using node_t = std::vector<node_data_t>;


using x3::float_;
using x3::eps;
using x3::lit;
using x3::_val;
using x3::_attr;    
using x3::attr;
using x3::alnum;
using ascii::char_;
using boost::fusion::at_c;

x3::rule<class label, std::string> const label = "label";
x3::rule<class ilabel, std::string> const ilabel = "ilabel";
x3::rule<class length, float> const length = "length";
x3::rule<class tip, node_t> const tip = "tip";
x3::rule<class inode, node_t> const inode = "inode";
x3::rule<class node, node_t> const node = "node";
x3::rule<class tree, node_t> const tree = "tree";

auto make_tip = [](auto& ctx) {
    auto s = at_c<0>(_attr(ctx));
    auto f = at_c<1>(_attr(ctx));
    _val(ctx) = node_t(1, {s,f});
};
auto make_inode = [](auto& ctx) {
    auto s = at_c<1>(_attr(ctx));
    auto f = at_c<2>(_attr(ctx));
    _val(ctx) = node_t(1, {s,f});
    _val(ctx).front().max_length = 0.0;
    auto v = at_c<0>(_attr(ctx));
    for(auto &&w : v) {
        auto n = _val(ctx).size();
        _val(ctx).front().total_length += w.front().total_length;
        _val(ctx).front().max_length = std::max(
            _val(ctx).front().max_length,
            f + w.front().max_length);
        for( auto &&x : w) {
            _val(ctx).push_back(x);
            _val(ctx).back().parent += n;
        }
        _val(ctx)[n].parent = 0;
    }
};

auto const tree_def = node >> -lit(';');
auto const node_def = tip | inode;
auto const tip_def = (label >> length)[make_tip];
auto const inode_def = (('(' >> (node % ',') >> ')') >> ilabel >> length)[make_inode];

auto const label_def = +(alnum | char_("-_%/."));
auto const ilabel_def = label | attr("");
auto const length_def = (':' >> float_) | attr(1.0);

BOOST_SPIRIT_DEFINE(label);
BOOST_SPIRIT_DEFINE(ilabel);
BOOST_SPIRIT_DEFINE(length);
BOOST_SPIRIT_DEFINE(tip);
BOOST_SPIRIT_DEFINE(inode);
BOOST_SPIRIT_DEFINE(node);
BOOST_SPIRIT_DEFINE(tree);

} // namespace newick

namespace mutk {
namespace detail { 

bool parse_newick(const std::string &text, pedigree_graph::Graph &graph,
    pedigree_graph::vertex_t root, bool normalize) {

    auto it = text.begin();
    auto end = text.end();

    newick::node_t phy;

    bool result = parse(it, end, newick::tree, phy);

    if(!(result && it == end)) {
        return false;
    }

    auto sex = get(boost::vertex_sex, graph, root);
    auto ploidy = get(boost::vertex_ploidy, graph, root);

    float scale = 1.0f;
    if(normalize) {
        scale = 1.0/phy[0].max_length;
    }

    std::vector<pedigree_graph::vertex_t> vertex_map(phy.size());
    for(std::size_t x = 0; x < phy.size(); ++x) {
        vertex_map[x] = add_vertex({phy[x].label, {sex, {ploidy, VertexType::Somatic}}}, graph);
        if(x == 0) {
            add_edge(root, vertex_map[0], {scale*phy[x].length, SOMA_EDGE}, graph);
        } else {
            add_edge(vertex_map[phy[x].parent], vertex_map[x], {scale*phy[x].length, SOMA_EDGE}, graph);
        }
    }

    return true;
}

TEST_CASE("[libmutk] detail::parse_newick") {
    using mutk::detail::parse_newick;
    using boost::edge_length;
    using boost::vertex_label;
    using boost::vertex_sex;
    using boost::vertex_ploidy;

    mutk::detail::pedigree_graph::Graph G;
    add_vertex({"Root", {Sex::Male, 2}},G);

    REQUIRE(parse_newick("(A:0.1,B:0.2)C_test:0.4", G, 0, false));
    REQUIRE(num_vertices(G) == 4);
    REQUIRE(num_edges(G) == 3);

    CHECK(get(vertex_label, G, 0) == "Root");
    CHECK(get(vertex_label, G, 1) == "C_test");
    CHECK(get(vertex_label, G, 2) == "A");
    CHECK(get(vertex_label, G, 3) == "B");
    CHECK(get(vertex_sex, G, 0) == Sex::Male);
    CHECK(get(vertex_sex, G, 1) == Sex::Male);
    CHECK(get(vertex_sex, G, 2) == Sex::Male);
    CHECK(get(vertex_sex, G, 3) == Sex::Male);
    CHECK(get(vertex_ploidy, G, 0) == 2);
    CHECK(get(vertex_ploidy, G, 1) == 2);
    CHECK(get(vertex_ploidy, G, 2) == 2);
    CHECK(get(vertex_ploidy, G, 3) == 2);

    REQUIRE(edge(0,1,G).second);
    REQUIRE(edge(1,2,G).second);
    REQUIRE(edge(1,3,G).second);

    CHECK(get(edge_length, G, edge(0,1,G).first) == 0.4f);
    CHECK(get(edge_length, G, edge(1,2,G).first) == 0.1f);
    CHECK(get(edge_length, G, edge(1,3,G).first) == 0.2f);

    REQUIRE(parse_newick("((S1,S2),S3)", G, 0, true));
    REQUIRE(num_vertices(G) == 9);
    REQUIRE(num_edges(G) == 8);

    REQUIRE(edge(0,4,G).second);
    REQUIRE(edge(4,5,G).second);
    REQUIRE(edge(4,8,G).second);
    REQUIRE(edge(5,6,G).second);
    REQUIRE(edge(5,7,G).second);

    CHECK(get(edge_length, G, edge(0,4,G).first) == 1.0f/3.0f);
    CHECK(get(edge_length, G, edge(4,5,G).first) == 1.0f/3.0f);
    CHECK(get(edge_length, G, edge(4,8,G).first) == 1.0f/3.0f);
    CHECK(get(edge_length, G, edge(5,6,G).first) == 1.0f/3.0f);
    CHECK(get(edge_length, G, edge(5,7,G).first) == 1.0f/3.0f);
}

}
}
