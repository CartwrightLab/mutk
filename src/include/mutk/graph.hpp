/*
# Copyright (c) 2023 Reed A. Cartwright <racartwright@gmail.com>
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

#ifndef MUTK_GRAPH_HPP
#define MUTK_GRAPH_HPP

#include "message.hpp"

#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>

// Install boost graph properties
namespace boost {
enum edge_length_t { edge_length };
enum edge_type_t { edge_type };

enum vertex_label_t { vertex_label };
enum vertex_data_t { vertex_data };
enum vertex_ploidy_t { vertex_ploidy };

BOOST_INSTALL_PROPERTY(edge, length);
BOOST_INSTALL_PROPERTY(edge, type);

BOOST_INSTALL_PROPERTY(vertex, label);
BOOST_INSTALL_PROPERTY(vertex, data);
BOOST_INSTALL_PROPERTY(vertex, ploidy);
}


namespace mutk {

// A strongly-type int. Use unitary + to do a static cast.
enum struct member_id_t : int {};
constexpr auto operator+(member_id_t value) {
    return static_cast<std::underlying_type_t<member_id_t>>(value);
}
enum struct sample_id_t : int {};
constexpr auto operator+(sample_id_t value) {
    return static_cast<std::underlying_type_t<sample_id_t>>(value);
}

namespace relationship_graph {

using EdgeLengthProp = boost::property<boost::edge_length_t, float>;
using EdgeProp = EdgeLengthProp;

using VertexDataProp   = boost::property<boost::vertex_data_t, std::vector<sample_id_t>>;
using VertexPloidyProp = boost::property<boost::vertex_ploidy_t, Ploidy, VertexDataProp>;
using VertexLabelProp = boost::property<boost::vertex_label_t, std::string, VertexPloidyProp>;
using VertexProp = VertexLabelProp;

} // namespace relationship_graph

namespace junction_tree {

using VertexLabelProp = boost::property<boost::vertex_label_t, std::vector<variable_t>>;
using VertexProp = VertexLabelProp;

} // namespace junctinon_tree

using RelationshipGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
        relationship_graph::VertexProp, relationship_graph::EdgeProp>;

using JunctionTree = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
        junction_tree::VertexProp, boost::no_property>;

template<class G>
auto make_vertex_range(G &graph) {
    return boost::make_iterator_range(vertices(graph));
}

template<class G>
auto make_adj_vertex_range(typename G::vertex_descriptor v,  G &graph) {
    return boost::make_iterator_range(adjacent_vertices(v, graph));
}

template<class G>
auto make_inv_vertex_range(typename G::vertex_descriptor v,  G &graph) {
    return boost::make_iterator_range(inv_adjacent_vertices(v, graph));
}

} // namespace mutk

#endif // MUTK_GRAPH_HPP
