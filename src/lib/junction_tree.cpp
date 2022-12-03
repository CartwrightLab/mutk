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

#include "unit_testing.hpp"

#include <boost/container/small_vector.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

#include "mutk/graph_builder.hpp"

using mutk::potential_t;
using mutk::clique_t;

using mutk::relationship_graph::JunctionTree;

using vertex_t = JunctionTree::vertex_descriptor;

using vertex_set_t = boost::container::flat_set<int, std::less<int>,
    boost::container::small_vector<int, 4>>;

template<class SinglePassRange1, class SinglePassRange2>
auto
set_union_size(const SinglePassRange1& rng1,
        const SinglePassRange2& rng2) {
    using size_type = typename std::iterator_traits<typename SinglePassRange1::const_iterator>::difference_type;
    size_type count = 0;
    auto it1 = std::begin(rng1);
    auto it2 = std::begin(rng2);
    for(; it1 != std::end(rng1) ; ++count ) {
        if(it2 == std::end(rng2)) {
            return count + std::distance(it1, std::end(rng1));
        }
        if(*it2 < *it1) {
            ++it2;
        } else {
            if(!(*it1 < *it2)) {
                ++it2;
            }
            ++it1;
        }
    }
    return count + std::distance(it2, std::end(rng2));
}

template<class SinglePassRange1, class SinglePassRange2>
auto
set_diff_size(const SinglePassRange1& rng1,
        const SinglePassRange2& rng2) {
    using size_type = typename std::iterator_traits<typename SinglePassRange1::const_iterator>::difference_type;
    size_type count = 0;
    auto it1 = std::begin(rng1);
    auto it2 = std::begin(rng2);
    while(it1 != std::end(rng1)) {
        if(it2 == std::end(rng2)) {
            return count + std::distance(it1, std::end(rng1));
        }
        if(*it1 < *it2) {
            ++count;
            ++it1;
        } else {
            if(!(*it2 < *it1)) {
                ++it1;
            }
            ++it2;
        }
    }
    return count;
}

namespace {
struct cost_result_t {
    size_t cost;
    bool   is_new_node;
    int    out_degree;
    const vertex_set_t * labels;
    vertex_t v;

    friend bool operator<(const cost_result_t & lhs, const cost_result_t & rhs) {
        return std::tie(lhs.cost, lhs.is_new_node, lhs.out_degree, *lhs.labels, lhs.v) <
               std::tie(rhs.cost, rhs.is_new_node, rhs.out_degree, *rhs.labels, rhs.v);
    }
};

class cost_visitor : public boost::default_dfs_visitor {
 public:
    cost_visitor(vertex_t query, cost_result_t &output, const std::vector<vertex_set_t> &node_labels) : 
        query_{query}, output_{output}, node_labels_(node_labels) {}

    void start_vertex(vertex_t v,  const JunctionTree & graph) {
        // initialize cost based on the union of v and query
        output_.cost = set_union_size(node_labels_[v], node_labels_[query_]);
        output_.is_new_node = true;
        output_.out_degree = out_degree(v, graph);
        output_.labels = &node_labels_[v];
        output_.v = v;
    }
    void finish_vertex(vertex_t v, const JunctionTree & graph) {
        // Unless we write our own DFS algorithm, includes must be called
        // twice per visited vertex. Once to check for termination.
        // Another to check if we have terminated.
        if(!boost::range::includes(node_labels_[v], node_labels_[query_])) {
            return;
        }
        // Prefer a vertex with the lowest cost
        size_t cost = node_labels_[v].size();
        if(cost > output_.cost) {
            return;
        }
        cost_result_t lhs;
        lhs.cost = cost;
        lhs.is_new_node = false;
        lhs.out_degree = out_degree(v, graph);
        lhs.labels = &node_labels_[v];
        lhs.v = v;
        if(lhs < output_) {
            output_ = lhs;
        }
    }

    vertex_t query_;
    cost_result_t &output_;
    const std::vector<vertex_set_t> &node_labels_;
};
}


// Try to link v2 into the graph that begins at v1
// If v2 is not a subset of v1, we will need to create a new node
// to link them. In this case the cost of the connection is the
// size of the union of v1 and v2, and the second element is empty.
// Otherwise, search v1 and all descendants looking for the smallest
// node that has v2 as a subset. If multiple nodes are equally small,
// prefer the one that has the smallest degree. If nodes are still tied,
// prefer the one with the smallest label set.
static
std::tuple<cost_result_t, vertex_t, vertex_t>
find_best_connection(vertex_t v1, vertex_t v2, const JunctionTree & graph,
        const std::vector<vertex_set_t> & node_labels) {
    if(node_labels[v1].size() < node_labels[v2].size()) {
        std::swap(v1,v2);
    }
    using Color = boost::color_traits< boost::default_color_type >;

    cost_result_t result;
    cost_visitor vis(v2, result, node_labels);

    std::vector< boost::default_color_type > color_map(
        num_vertices(graph), Color::white());

    boost::default_color_type c = Color::white();

    auto term_func = [&](vertex_t v, const JunctionTree & g) {
        return !boost::range::includes(node_labels[v], node_labels[v2]);
    };

    boost::depth_first_visit(graph, v1, vis,
        boost::make_iterator_property_map(color_map.begin(),
            get(boost::vertex_index, graph), c), term_func);

    return {result, v1, v2};
}

mutk::relationship_graph::JunctionTree
create_junction_tree(const mutk::relationship_graph::Graph &graph,
    const std::vector<potential_t> &potentials,
    const std::vector<clique_t> &elimination_order) {

    using mutk::relationship_graph::Graph;

    using mutk::relationship_graph::variable_t;

    // a 'variable' is a vertex in the relationship graph.
    // a 'rank' is the position of the variable in the elimination order
    // variables eliminated first have lower rank
    std::vector<variable_t> rank_to_var(num_vertices(graph));
    std::vector<int> var_to_rank(num_vertices(graph));
    for(int i=0;i<elimination_order.size();++i) {
        Graph::vertex_descriptor v = elimination_order[i].front();
        rank_to_var[i] = variable_t(v);
        var_to_rank[v] = i;
    }

    // Each node in a junction tree is labeled with sets of variables.
    // While constructing the junction tree we will use ranks and
    // external storage.
    mutk::relationship_graph::JunctionTree junction(potentials.size());

    // Initialize nodes based on potentials.
    // Set all nodes to active
    std::vector<vertex_set_t> node_labels(potentials.size());
    std::vector<vertex_t> active_nodes, working_nodes;

    for(int i=0; i < potentials.size(); ++i) {
        for(auto x : potentials[i].variables ) {
            node_labels[i].insert(var_to_rank[+x]);
        }
        active_nodes.push_back(i);
    }

    // Construct the augmented junction tree using binary fusion
    for(int var=0; var < elimination_order.size(); ++var) {
        // identify active nodes that contain `var` and move them to working nodes
        working_nodes.clear();
        boost::remove_erase_if(active_nodes, [&](auto v){
            if(node_labels[v].contains(var)) {
                working_nodes.push_back(v);
                return true;
            }
            return false;
        });
        while(working_nodes.size() > 1) {
            // find which pair of labels can be merged with lowest cost
            auto it = working_nodes.begin();
            auto jt = std::next(it);
            auto bit = it;
            auto bjt = jt;
            auto best = find_best_connection(*it, *jt, junction, node_labels);
            for(++jt; jt != working_nodes.end(); ++it) {
                for(it = working_nodes.begin(); it != jt; ++it) {
                    auto result = find_best_connection(*it, *jt, junction, node_labels);
                    if(std::get<0>(result) < std::get<0>(best)) {
                        best = result;
                        bjt = jt;
                        bit = it;
                    }
                }
            }
            // Result holds information about the best location to attach v2 to v1.
            // bit and bjt are iterators to the locations of v2 and v1. bit occurs
            // before bjt in the sequence. bit and bjt can refer to either v1 or v2.
            auto [result, v1, v2] = best;
            if(result.is_new_node) {
                // construct a new node and label it with the
                // union of the two nodes.
                auto u = add_vertex(junction);
                node_labels.push_back(node_labels[v1]);
                node_labels.back().merge(node_labels[v2]);
                add_edge(u,v1,junction);
                add_edge(u,v2,junction);
                *bit = u;
            } else {
                // attach v2 to u, duplicating node as needed
                auto u = result.v;
                if(out_degree(u, junction) > 1) {
                    auto w = add_vertex(junction);
                    node_labels.push_back(node_labels[u]);
                    std::vector<vertex_t> adj;
                    auto [vi, vi_end] = adjacent_vertices(u, junction);
                    for(; vi != vi_end; ++vi) {
                        adj.push_back(*vi);
                    }
                    for(auto v : adj) {
                        add_edge(w, v, junction);
                    }
                    clear_out_edges(u, junction);
                    add_edge(u, w, junction);

                }
                add_edge(u, v2, junction);
                *bit = v1;
            }
            // bit refers to a new node
            // bjt needs to be erased
            working_nodes.erase(bjt);
        }
        // only 1 working node left
        vertex_t v = working_nodes.front();
        if(node_labels[v].size() > 1) {
            vertex_t u = add_vertex(junction);
            add_edge(u, v, junction);
            // eliminate var from the label of the new node
            node_labels.push_back(node_labels[v]);
            node_labels.back().erase(var);
            active_nodes.push_back(u);
        }
    }
    return junction;
}