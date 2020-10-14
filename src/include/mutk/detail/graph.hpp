/*
# Copyright (c) 2014-2020 Reed A. Cartwright <reed@cartwright.ht>
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

#ifndef MUTK_DETAIL_GRAPH_HPP
#define MUTK_DETAIL_GRAPH_HPP

#include <mutk/pedigree.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/container/flat_set.hpp>

// Install boost graph properties
namespace boost {
enum edge_length_t { edge_length };
enum edge_type_t { edge_type };

enum vertex_sex_t { vertex_sex };
enum vertex_ploidy_t {vertex_ploidy };
enum vertex_label_t { vertex_label };
enum vertex_type_t { vertex_type };

enum edge_family_t { edge_family };
enum vertex_group_t { vertex_group };

BOOST_INSTALL_PROPERTY(edge, length);
BOOST_INSTALL_PROPERTY(edge, type);
BOOST_INSTALL_PROPERTY(vertex, sex);
BOOST_INSTALL_PROPERTY(vertex, ploidy);
BOOST_INSTALL_PROPERTY(vertex, label);
BOOST_INSTALL_PROPERTY(vertex, type);

BOOST_INSTALL_PROPERTY(edge, family);
BOOST_INSTALL_PROPERTY(vertex, group);
}

namespace mutk {
namespace detail {

enum EdgeType : int {
    GERM_EDGE = 0x1,
    SOMA_EDGE = 0x2,
};

enum struct VertexType : int {
    Germline, Somatic, Sample
};

using Sex = mutk::Pedigree::Sex;

namespace pedigree_graph {
// a member graph is a directed graph that we construct while parsing
// the pedigree data

using EdgeTypeProp = boost::property<boost::edge_type_t, EdgeType>;
using EdgeLengthProp = boost::property<boost::edge_length_t, float, EdgeTypeProp>;
using EdgeProp = EdgeLengthProp;

using VertexTypeProp = boost::property<boost::vertex_type_t, VertexType>;
using VertexPloidyProp = boost::property<boost::vertex_ploidy_t, int, VertexTypeProp>;
using VertexSexProp = boost::property<boost::vertex_sex_t, Sex, VertexPloidyProp>;
using VertexLabelProp = boost::property<boost::vertex_label_t, std::string, VertexSexProp>;
using VertexProp = VertexLabelProp;

using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
        VertexProp, EdgeProp>;
using vertex_t = boost::graph_traits<Graph>::vertex_descriptor;
using edge_t = boost::graph_traits<Graph>::edge_descriptor;

static_assert(std::is_integral<vertex_t>::value,
    "vertex_t is not an integral type, this violates many assumptions that have been made.");

inline void print_graph(const Graph &graph, std::ostream& os = std::cout) {
    for(auto [vi, ve] = vertices(graph); vi != ve; ++vi) {
        os << *vi
           << "\t" << get(boost::vertex_label, graph, *vi)
           << "\t" << get(boost::vertex_ploidy, graph, *vi)
           << "\n";
    }
    for(auto [ei, ee] = edges(graph); ei != ee; ++ei) {
        os << source(*ei,graph) << "-->" << target(*ei,graph) << "\n";
    }
}

}

bool parse_newick(const std::string &text,
    pedigree_graph::Graph &graph,
    pedigree_graph::vertex_t root,
    bool normalize);

namespace junction_tree {
template<typename V_, int N_>
using small_vector_t = boost::container::small_vector<V_, N_>;
//using small_vector_t = std::vector<pedigree_graph::vertex_t>;

template<typename V_, int N_>
using flat_set_t = boost::container::flat_set<V_, std::less<V_>, small_vector_t<V_, N_>>;
using neighbors_t = flat_set_t<pedigree_graph::vertex_t, 4>;

using VertexProp = boost::property<boost::vertex_label_t, neighbors_t>;

using EdgeProp = boost::property<boost::edge_color_t, boost::default_color_type>;

using Graph = boost::adjacency_list<boost::vecS, boost::vecS,
    boost::undirectedS, VertexProp, EdgeProp>;

using vertex_t = boost::graph_traits<Graph>::vertex_descriptor;
using edge_t = boost::graph_traits<Graph>::edge_descriptor;
}



} // namespace mutk::detail
} // namespace mutk

#endif // MUTK_DETAIL_GRAPH_HPP
