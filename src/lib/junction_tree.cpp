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

template<class G>
auto make_vertex_range(G &graph) {
    return boost::make_iterator_range(vertices(graph));
}

using mutk::potential_t;
using mutk::clique_t;

using mutk::relationship_graph::JunctionTree;

using vertex_t = JunctionTree::vertex_descriptor;

using vertex_set_t = boost::container::flat_set<int, std::less<int>,
    boost::container::small_vector<int, 4>>;

namespace {
struct source_t {
    vertex_t source;
    int label_size;
    bool not_a_subset_of_source;
    int out_degree;
    vertex_set_t label;

    friend bool operator<(const source_t & lhs, const source_t & rhs) {
        return std::make_tuple(lhs.label_size, lhs.not_a_subset_of_source, lhs.out_degree, lhs.label) <
               std::make_tuple(rhs.label_size, rhs.not_a_subset_of_source, rhs.out_degree, rhs.label);
    }
};
} // anon namespace

// Try to link v2 into the graph that begins at v1
// If v2 is not a subset of v1, we will need to create a new node
// to link them. In this case the cost of the connection is the
// size of the union of v1 and v2, and the second element is empty.
// Otherwise, search v1 and all descendants looking for the smallest
// node that has v2 as a subset. If multiple nodes are equally small,
// prefer the one that has the smallest degree. If nodes are still tied,
// prefer the one with the smallest label set.
static
source_t
find_best_connection(vertex_t root, vertex_t query, const JunctionTree & graph,
        const std::vector<vertex_set_t> & node_labels) {
    using boost::range::includes;

    // attach query to a node that is the union of root and query
    source_t result;
    result.source = root;
    result.label = node_labels[root];
    result.label.insert(boost::container::ordered_unique_range_t{},
        node_labels[query].begin(), node_labels[query].end());
    result.label_size = result.label.size();
    result.not_a_subset_of_source = true;
    result.out_degree = 1;

    // If query is not a subset of root node then the only place we can
    // attach is at a new node.
    if(!includes(node_labels[root], node_labels[query])) {
        return result;
    }
    // do a depth first search looking for a better attachment point
    using AdjIter = JunctionTree::adjacency_iterator;
    using VertexInfo = std::pair<vertex_t, std::pair<AdjIter, AdjIter>>;
    std::vector<VertexInfo> stack;

    // initialize our search at the root node
    vertex_t u = root;
    auto [ai, ai_end] = adjacent_vertices(u, graph);
    stack.push_back({u, {ai, ai_end}});
    while(!stack.empty()) {
        u = stack.back().first;
        std::tie(ai, ai_end) = stack.back().second;
        stack.pop_back();
        while(ai != ai_end) {
            // Change context of search if query is a subset of an adjacent vertex
            vertex_t v = *ai;
            if(includes(node_labels[v], node_labels[query])) {
                stack.push_back({u, {++ai, ai_end}});
                u = v;
                std::tie(ai, ai_end) = adjacent_vertices(u, graph);
            } else {
                ++ai;
            }
        }
        // Test if this node is a better match than the current best node.
        if(node_labels[u].size() <= result.label_size) {
            source_t here;
            here.source = u;
            here.label = node_labels[u];
            here.label_size = here.label.size();
            here.not_a_subset_of_source = false;
            here.out_degree = out_degree(u, graph);
            if(here < result) {
                result = here;
            }
        }
    }
    return result;
}

mutk::relationship_graph::JunctionTree
mutk::create_junction_tree(const mutk::relationship_graph::Graph &graph,
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

    // helper lambda function
    auto order_nodes = [&](vertex_t v1, vertex_t v2) {
        if(node_labels[v1].size() < node_labels[v2].size()) {
            std::swap(v1,v2);
        }
        return std::make_pair(v1,v2);
    };

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
            auto [v1, v2] = order_nodes(*it, *jt);
            auto best = find_best_connection(v1, v2, junction, node_labels);
            auto best_pair = std::make_pair(it, jt);
            for(++jt; jt != working_nodes.end(); ++it) {
                for(it = working_nodes.begin(); it != jt; ++it) {
                    auto [v1, v2] = order_nodes(*it, *jt);
                    auto result = find_best_connection(v1, v2, junction, node_labels);
                    if(result < best) {
                        best = result;
                        best_pair = {it, jt};
                    }
                }
            }
            // Result holds information about the best location to attach v2 to v1.
            // bit and bjt are iterators to the locations of v2 and v1. bit occurs
            // before bjt in the sequence. bit and bjt can refer to either v1 or v2.
            std::tie(it, jt) = best_pair;
            std::tie(v1, v2) = order_nodes(*it, *jt);
            if(best.not_a_subset_of_source) {
                // construct a new node and label it with the
                // union of the two nodes.
                auto u = add_vertex(junction);
                node_labels.push_back(node_labels[v1]);
                node_labels.back().insert(boost::container::ordered_unique_range_t{},
                    node_labels[v2].begin(), node_labels[v2].end());
                add_edge(u,v1,junction);
                add_edge(u,v2,junction);
                *it = u;
            } else {
                // attach v2 to u, duplicating node as needed
                auto u = best.source;
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
                *it = v1;
            }
            // it refers to a new node
            // jt needs to be erased
            working_nodes.erase(jt);
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

    // Label Junction Tree
    for(auto v : make_vertex_range(junction)) {
        std::vector<variable_t> vars;
        for(auto rank : node_labels[v]) {
            vars.push_back(rank_to_var[rank]);
        }
        put(boost::vertex_label, junction, v, vars);
    }

    return junction;
}

// LCOV_EXCL_START
TEST_CASE("create_junction_tree() constructs a junction tree.") {
    using mutk::relationship_graph::Graph;
    using mutk::relationship_graph::variable_t;

    auto mkpot = [](const std::initializer_list<int> & range) -> potential_t {
        std::vector<variable_t> ret;
        for(auto v : range) {
            ret.push_back(variable_t(v));
        }
        return {ret, {}};
    };

    {
        Graph graph(7);
        add_edge(0,2,graph);
        add_edge(1,2,graph);
        add_edge(2,3,graph);
        add_edge(2,4,graph);
        add_edge(5,6,graph);

        std::vector<potential_t> potentials;
        potentials.push_back(mkpot({0}));
        potentials.push_back(mkpot({1}));
        potentials.push_back(mkpot({2}));
        potentials.push_back(mkpot({2,0,1}));
        potentials.push_back(mkpot({3}));
        potentials.push_back(mkpot({3,2}));
        potentials.push_back(mkpot({4}));
        potentials.push_back(mkpot({4,2}));
        potentials.push_back(mkpot({5}));
        potentials.push_back(mkpot({6}));
        potentials.push_back(mkpot({6,5}));

        std::vector<clique_t> cliques;
        cliques.push_back(clique_t({6,5}));
        cliques.push_back(clique_t({5}));
        cliques.push_back(clique_t({4,2}));
        cliques.push_back(clique_t({3,2}));
        cliques.push_back(clique_t({2,0,1}));
        cliques.push_back(clique_t({1,0}));
        cliques.push_back(clique_t({0}));

        // ex 6: 6 <- 6,5 ; 6,5 <- 5
        // ex 5: 
        // ex 4: 4 <- 4,2 ; 4,2 <- 2
        // ex 3: 3 <- 3,2 ; 3,2 <- 3
        // ex 2: 2 <- 2,0,1 ; 2,0,1 <- 0,1
        // ex 1: 0,1 <- 0 
        // ex 0:
        auto tree = create_junction_tree(graph, potentials, cliques);

        CHECK(num_vertices(tree) == 12);

    }
}
// LCOVE_EXCL_STOP