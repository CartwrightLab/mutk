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

#include <iostream>

#include <boost/container/small_vector.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/graph/topological_sort.hpp>

#include "mutk/graph_builder.hpp"

template<class G>
auto make_vertex_range(G &graph) {
    return boost::make_iterator_range(vertices(graph));
}

template<class G>
auto make_adj_vertex_range(typename G::vertex_descriptor v,  G &graph) {
    return boost::make_iterator_range(adjacent_vertices(v, graph));
}

using mutk::xpotential_t;
using mutk::clique_t;



using vertex_set_t = boost::container::flat_set<int, std::less<int>,
    boost::container::small_vector<int, 4>>;

static
vertex_set_t merge_vertex_sets(vertex_set_t a, const vertex_set_t &b) {
    a.insert(boost::container::ordered_unique_range_t{},
        b.begin(), b.end());
    return a;
}

namespace {
using VertexLabelProp = boost::property<boost::vertex_label_t, vertex_set_t>;
using VertexProp = VertexLabelProp;

using directed_tree_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
        VertexProp, boost::no_property>;

using vertex_t = directed_tree_t::vertex_descriptor;

struct source_t {
    vertex_t source;
    int label_size;
    bool creates_new_node;
    int out_degree;
    vertex_set_t label;

    friend bool operator<(const source_t & lhs, const source_t & rhs) {
        return std::make_tuple(lhs.label_size, lhs.creates_new_node, lhs.out_degree, lhs.label) <
               std::make_tuple(rhs.label_size, rhs.creates_new_node, rhs.out_degree, rhs.label);
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
find_best_connection(vertex_t root, vertex_t query, const directed_tree_t & graph) {
    using boost::range::includes;

    auto labels = get(boost::vertex_label, graph);

    // attach query to a node that is the union of root and query
    source_t result;
    result.source = root; 
    result.label = merge_vertex_sets(labels[root], labels[query]);
    result.label_size = result.label.size();
    result.creates_new_node = true;
    result.out_degree = 1;

    // If query is not a subset of root node then the only place we can
    // attach is at a new node.
    if(!includes(labels[root], labels[query])) {
        return result;
    }
    // do a depth first search looking for a better attachment point
    using AdjIter = directed_tree_t::adjacency_iterator;
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
            if(includes(labels[v], labels[query])) {
                stack.push_back({u, {++ai, ai_end}});
                u = v;
                std::tie(ai, ai_end) = adjacent_vertices(u, graph);
            } else {
                ++ai;
            }
        }
        // Test if this node is a better match than the current best node.
        if((int)labels[u].size() <= result.label_size) {
            source_t here;
            here.source = u;
            here.label = labels[u];
            here.label_size = here.label.size();
            here.creates_new_node = false;
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
    const std::vector<xpotential_t> &potentials,
    const std::vector<clique_t> &elimination_order) {

    using mutk::relationship_graph::variable_t;

    // a 'variable' is a vertex in the relationship graph.
    // a 'rank' is the position of the variable in the elimination order
    // variables eliminated first have lower rank
    std::vector<variable_t> rank_to_var(num_vertices(graph));
    std::vector<int> var_to_rank(num_vertices(graph));
    for(size_t i = 0; i < elimination_order.size(); ++i) {
        auto v = elimination_order[i].front();
        rank_to_var[i] = variable_t(v);
        var_to_rank[v] = i;
    }

    // Each node in a junction tree is labeled with sets of variables.
    // While constructing the junction tree we will use ranks and
    // external storage.
    directed_tree_t active_tree;

    // Initialize nodes based on potentials.
    // Set all nodes to active
    std::map<vertex_set_t, vertex_t> nodes;

    auto labels = get(boost::vertex_label, active_tree);

    std::vector<vertex_t> active_nodes, working_nodes;

    for(size_t i=0; i < potentials.size(); ++i) {
        vertex_set_t labels;
        for(auto x : potentials[i].variables ) {
            labels.insert(var_to_rank[+x]);
        }
        auto v = add_vertex(labels, active_tree);
        nodes.try_emplace(labels, v);
        active_nodes.push_back(v);
    }

    // helper lambda function
    auto order_nodes = [&](vertex_t v1, vertex_t v2) {
        if(labels[v1].size() < labels[v2].size()) {
            std::swap(v1,v2);
        }
        return std::make_pair(v1,v2);
    };

    // Construct the augmented junction tree using binary fusion
    for(size_t var=0; var < elimination_order.size(); ++var) {
        // identify active nodes that contain `var` and move them to working nodes
        working_nodes.clear();
        boost::remove_erase_if(active_nodes, [&](auto v){
            if(labels[v].contains(var)) {
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
            auto best = find_best_connection(v1, v2, active_tree);
            auto best_pair = std::make_pair(it, jt);
            for(++jt; jt != working_nodes.end(); ++jt) {
                for(it = working_nodes.begin(); it != jt; ++it) {
                    auto [v1, v2] = order_nodes(*it, *jt);
                    auto result = find_best_connection(v1, v2, active_tree);
                    if(result < best) {
                        best = result;
                        best_pair = {it, jt};
                    }
                }
            }
            // `best` holds information about the best location to attach v2 to v1.
            // bit and bjt are iterators to the locations of v2 and v1. bit occurs
            // before bjt in the sequence. bit and bjt can refer to either v1 or v2.
            std::tie(it, jt) = best_pair;
            std::tie(v1, v2) = order_nodes(*it, *jt);

            if(best.creates_new_node) {
                vertex_t u;
                // check if node already exists and create it if necessary
                auto node_it = nodes.find(best.label);
                if(node_it == nodes.end()) {
                    u = add_vertex(best.label, active_tree);
                    nodes.try_emplace(best.label, u);
                } else {
                    u = node_it->second;
                }
                add_edge(u, v1, active_tree);
                add_edge(u, v2, active_tree);
                *it = u;
            } else {
                // attach v2 to u
                auto u = best.source;
                add_edge(u, v2, active_tree);
                *it = v1;
            }
            // it refers to a new node
            // jt needs to be erased
            working_nodes.erase(jt);
        }
        // only 1 working node left
        vertex_t v = working_nodes.front();
        if(labels[v].size() > 1) {
            // Create a new node by eliminating var from v's label
            auto label = labels[v];
            label.erase(var);
            // check if node already exists and create it if necessary
            auto node_it = nodes.find(label);
            vertex_t u;
            if(node_it == nodes.end()) {
                u = add_vertex(label, active_tree);
                nodes.try_emplace(label, u);
                // push the new node onto active_nodes
                active_nodes.push_back(u);
            } else {
                u = node_it->second;
            } 
            add_edge(u, v, active_tree);
        }
    }

    // Create Binary Tree
    for(auto v : make_vertex_range(active_tree)) {
        if(out_degree(v, active_tree) < 3) {
            continue;
        }

        // Remember which vertices are connected to v.
        auto adj_range = make_adj_vertex_range(v, active_tree);
        std::vector<vertex_t> adj_vertices;
        boost::copy(adj_range, std::back_inserter(adj_vertices));
        clear_out_edges(v, active_tree);

        // construct a caterpillar tree.
        auto u = v;
        size_t i;
        for(i=0; i < adj_vertices.size()-2;++i) {
            auto w = add_vertex(labels[v], active_tree);
            add_edge(u, w, active_tree);
            add_edge(u, adj_vertices[i], active_tree);
            u = w;
        }
        add_edge(u, adj_vertices[i], active_tree);
        add_edge(u, adj_vertices[i+1], active_tree);
    }

    // Sort by topological order
    std::vector<vertex_t> rev_topo_order;
    topological_sort(active_tree, std::back_inserter(rev_topo_order));

    std::vector<int> vertex_to_order(rev_topo_order.size());
    for(size_t i = 0; i < rev_topo_order.size(); ++i) {
        vertex_to_order[rev_topo_order[i]] = i;
    }

    // Create Junction Tree in reverse topological order
    mutk::relationship_graph::JunctionTree output(num_vertices(active_tree));
    // Create Edges
    for(auto [ei, ei_end] = edges(active_tree); ei != ei_end; ++ei) {
        auto s = source(*ei, active_tree);
        auto t = target(*ei, active_tree);
        add_edge(vertex_to_order[s], vertex_to_order[t], output);
    }

    // Label Junction Tree Vertices
    for(auto v : rev_topo_order) {
        std::vector<variable_t> vars;
        for(auto i : labels[v]) {
            vars.push_back(rank_to_var[i]);
        }
        put(boost::vertex_label, output, vertex_to_order[v], vars);
    }

    return output;
}

// LCOV_EXCL_START
TEST_CASE("create_junction_tree() constructs a junction tree.") {
    using mutk::relationship_graph::Graph;
    using mutk::relationship_graph::JunctionTree;
    using mutk::relationship_graph::variable_t;

    auto mkpot = [](const std::initializer_list<int> & range) -> xpotential_t {
        std::vector<variable_t> ret;
        for(auto v : range) {
            ret.push_back(variable_t(v));
        }
        return {ret, {}};
    };

    auto get_label = [](JunctionTree::vertex_descriptor a, const JunctionTree &g) {
        auto lab = get(boost::vertex_label, g, a);
        std::string output = std::to_string(+lab[0]);
        for(size_t i=1;i<lab.size();++i) {
            output += ",";
            output += std::to_string(+lab[i]);
        }
        return output;
    };

    auto get_labels = [&](const JunctionTree &g) {
        std::vector<std::string> ret;
        for(auto v : make_vertex_range(g)) {
            ret.push_back(get_label(v, g));
        }
        return ret;
    };

    auto get_edges = [&](const JunctionTree &g) {
        std::vector<std::string> ret;
        for(auto [ei, ei_end] = edges(g); ei != ei_end; ++ei) {
            std::string s = get_label(source(*ei, g), g);
            std::string t = get_label(target(*ei, g), g);
            s += "->";
            s += t;
            ret.push_back(std::move(s));
        }

        return ret;
    };

    {
        Graph graph(7);
        add_edge(0,2,graph);
        add_edge(1,2,graph);
        add_edge(2,3,graph);
        add_edge(2,4,graph);
        add_edge(5,6,graph);

        std::vector<xpotential_t> potentials;
        potentials.push_back(mkpot({2,0,1}));  // 0
        potentials.push_back(mkpot({3,2}));    // 1
        potentials.push_back(mkpot({4,2}));    // 2
        potentials.push_back(mkpot({6,5}));    // 3
        potentials.push_back(mkpot({0}));      // 4
        potentials.push_back(mkpot({1}));      // 5
        potentials.push_back(mkpot({2}));      // 6
        potentials.push_back(mkpot({3}));      // 7
        potentials.push_back(mkpot({4}));      // 8
        potentials.push_back(mkpot({5}));      // 9
        potentials.push_back(mkpot({6}));      //10

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

        std::vector<std::string> expected_labels = {
            "4", "4,2", "3", "3,2", "2", "2,1,0",
            "6", "6,5", "1", "1,0", "0", "5"};
        CHECK_EQ_RANGES(get_labels(tree), expected_labels);

        std::vector<std::string> expected_edges = {
            "4,2->4", "3,2->3", "2->4,2", "2->3,2", "2,1,0->2",
            "6,5->6", "1,0->2,1,0", "1,0->1", "0->1,0", "5->6,5"
        };
        CHECK_EQ_RANGES(get_edges(tree), expected_edges);
    }
    {
        Graph graph(5);
        add_edge(0,1,graph);
        add_edge(0,2,graph);
        add_edge(0,3,graph);
        add_edge(0,4,graph);

        std::vector<xpotential_t> potentials;
        potentials.push_back(mkpot({1,0}));    // 0
        potentials.push_back(mkpot({2,0}));    // 1
        potentials.push_back(mkpot({3,0}));    // 2
        potentials.push_back(mkpot({4,0}));    // 3
        potentials.push_back(mkpot({0}));      // 4
        potentials.push_back(mkpot({1}));      // 5
        potentials.push_back(mkpot({2}));      // 6
        potentials.push_back(mkpot({3}));      // 7
        potentials.push_back(mkpot({4}));      // 8

        std::vector<clique_t> cliques;
        cliques.push_back(clique_t({4,0}));
        cliques.push_back(clique_t({3,0}));
        cliques.push_back(clique_t({2,0}));
        cliques.push_back(clique_t({1,0}));
        cliques.push_back(clique_t({0}));

        auto tree = create_junction_tree(graph, potentials, cliques);

        std::vector<std::string> expected_labels = {
            "1", "1,0", "2", "2,0", "3", "3,0",
            "4", "4,0", "0", "0", "0"};
        CHECK_EQ_RANGES(get_labels(tree), expected_labels);

        std::vector<std::string> expected_edges = {
            "1,0->1", "2,0->2", "3,0->3", "4,0->4", "0->2,0",
            "0->1,0", "0->0", "0->3,0", "0->0", "0->4,0"
        };
        CHECK_EQ_RANGES(get_edges(tree), expected_edges);
    }
    {   // FIRST COUSIN MARRIAGE
        Graph graph(9);
        add_edge(0,3,graph);
        add_edge(1,3,graph);
        add_edge(0,4,graph);
        add_edge(1,4,graph);
        add_edge(2,6,graph);
        add_edge(3,6,graph);
        add_edge(4,7,graph);
        add_edge(5,7,graph);
        add_edge(6,8,graph);
        add_edge(7,8,graph);

        std::vector<xpotential_t> potentials;
        potentials.push_back(mkpot({3,0,1}));
        potentials.push_back(mkpot({4,0,1}));
        potentials.push_back(mkpot({6,2,3}));
        potentials.push_back(mkpot({7,4,5}));
        potentials.push_back(mkpot({8,6,7}));
        potentials.push_back(mkpot({0}));
        potentials.push_back(mkpot({1}));
        potentials.push_back(mkpot({2}));
        potentials.push_back(mkpot({3}));
        potentials.push_back(mkpot({4}));
        potentials.push_back(mkpot({5}));
        potentials.push_back(mkpot({6}));
        potentials.push_back(mkpot({7}));
        potentials.push_back(mkpot({8}));

        std::vector<clique_t> cliques;
        cliques.push_back(clique_t({8,6,7}));
        cliques.push_back(clique_t({5,4,7}));
        cliques.push_back(clique_t({2,3,6}));
        cliques.push_back(clique_t({7,4,6}));
        cliques.push_back(clique_t({6,3,4}));
        cliques.push_back(clique_t({4,0,1,3}));
        cliques.push_back(clique_t({3,0,1}));
        cliques.push_back(clique_t({1,0}));
        cliques.push_back(clique_t({0}));

        auto tree = create_junction_tree(graph, potentials, cliques);

        std::vector<std::string> expected_labels = {
            "4,1,0", "8", "8,7,6", "7", "7,6", "5", "5,7,4", "4",
            "7,4", "7,6,4", "6", "6,4", "2", "2,6,3", "3", "6,3",
            "6,4,3", "4,3", "4,3,1,0", "3,1,0", "1", "1,0", "0",
        };
        CHECK_EQ_RANGES(get_labels(tree), expected_labels);

        std::vector<std::string> expected_edges = {
            "8,7,6->8", "7,6->8,7,6", "7,6->7", "5,7,4->5", "7,4->5,7,4", "7,4->4",
            "7,6,4->7,6", "7,6,4->7,4", "6,4->7,6,4", "6,4->6", "2,6,3->2", "6,3->2,6,3",
            "6,3->3", "6,4,3->6,4", "6,4,3->6,3", "4,3->6,4,3", "4,3,1,0->4,1,0",
            "4,3,1,0->4,3", "3,1,0->4,3,1,0", "1,0->3,1,0", "1,0->1", "0->1,0"
        };
        CHECK_EQ_RANGES(get_edges(tree), expected_edges);        
    }
}
// LCOVE_EXCL_STOP
