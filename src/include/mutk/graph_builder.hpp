/*
# Copyright (c) 2022 Reed A. Cartwright <racartwright@gmail.com>
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

#ifndef MUTK_GRAPH_BUILDER_HPP
#define MUTK_GRAPH_BUILDER_HPP

#include <string>
#include <vector>
#include <unordered_map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>

#include "inheritance_model.hpp"

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

enum struct VertexPloidy : int {
    Noploid = 0,
    Haploid = 1,
    Diploid = 2
};

enum struct variable_t : int {};
constexpr auto operator+(variable_t value) {
    return static_cast<std::underlying_type_t<variable_t>>(value);
}

namespace graph {

using EdgeLengthProp = boost::property<boost::edge_length_t, float>;
using EdgeProp = EdgeLengthProp;

using VertexDataProp   = boost::property<boost::vertex_data_t, std::vector<sample_id_t>>;
using VertexPloidyProp = boost::property<boost::vertex_ploidy_t, VertexPloidy, VertexDataProp>;
using VertexLabelProp = boost::property<boost::vertex_label_t, std::string, VertexPloidyProp>;
using VertexProp = VertexLabelProp;

using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
        VertexProp, EdgeProp>;
}

namespace junction_tree {

using VertexLabelProp = boost::property<boost::vertex_label_t, std::vector<variable_t>>;
using VertexProp = VertexLabelProp;

using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
        VertexProp, boost::no_property>;
}

using Graph = graph::Graph;
using JunctionTree = junction_tree::Graph;

}

struct potential_t {
    std::vector<relationship_graph::variable_t> variables;
    std::vector<float> edge_lengths;
};

using clique_t = std::vector<mutk::relationship_graph::Graph::vertex_descriptor>;

class GraphBuilder {
 public:
    GraphBuilder();

    // Add families as biconnected components
    member_id_t AddSingle(const std::string &name, const std::string &sex, const std::vector<std::string> &sample_names);

    member_id_t AddPair(const std::string &name, const std::string &sex, const std::vector<std::string> &sample_names,
        const std::string &parent_name, float mutation_scale);

    member_id_t AddTrio(const std::string &name, const std::string &sex, const std::vector<std::string> &sample_names,
        const std::string &parent_name_a, float mutation_scale_a,
        const std::string &parent_name_b, float mutation_scale_b);

    // Add known samples
    void AddSamples(const std::vector<std::string> &sample_names);

    // Build the final relationship graph based on an inhertiance model
    void BuildGraph(const InheritanceModel &model, float mu);

 private:
    member_id_t LookupName(const std::string &name);

    relationship_graph::Graph CreateInitialGraph(const InheritanceModel &model, float mu);

    struct component_t {
        std::vector<member_id_t> members;
        std::vector<float> scales;
    };

    int member_counter_{0};

    std::unordered_map<std::string, member_id_t> map_member_name_to_id_;
    std::vector<std::string> member_names_;
    std::vector<std::string> member_sexes_;
    std::vector<std::vector<std::string>> member_input_samples_;
    std::vector<std::vector<sample_id_t>> member_data_samples_;

    std::unordered_map<std::string, sample_id_t> map_sample_name_to_id_;
    std::vector<std::string> sample_names_;

    std::vector<component_t> components_;
};


mutk::relationship_graph::JunctionTree
create_junction_tree(const mutk::relationship_graph::Graph &graph,
    const std::vector<potential_t> &potentials,
    const std::vector<clique_t> &elimination_order);

}

#endif // MUTK_GRAPH_BUILDER_HPP
