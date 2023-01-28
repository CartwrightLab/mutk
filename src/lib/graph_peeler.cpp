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
#include "unit_testing.hpp"


#include <mutk/graph.hpp>
#include <mutk/graph_peeler.hpp>

#include <boost/heap/d_ary_heap.hpp>

#include "junction_tree.hpp"

using mutk::junction_tree::component_t;
using mutk::junction_tree::clique_t;
using mutk::make_adj_vertex_range;
using mutk::make_vertex_range;
using mutk::make_inv_vertex_range;

static std::vector<clique_t>
triangulate_graph(const mutk::RelationshipGraph &graph);

static std::vector<component_t>
calculate_components(const mutk::RelationshipGraph &graph);

mutk::GraphPeeler mutk::GraphPeeler::Create(mutk::RelationshipGraph graph) {
    GraphPeeler peeler;

    peeler.graph_ = std::move(graph);

    auto components = calculate_components(peeler.graph_);
    auto cliques = triangulate_graph(peeler.graph_);

    peeler.tree_ = create_junction_tree(peeler.graph_, components, cliques);

    return peeler;
}

// Triangulate a graph where the vertices are in topological order
//
// Almond and Kong (1991) Optimality Issues in Constructing a Markov Tree from Graphical Models.
//     Research Report 329. University of Chicago, Dept. of Statistics
static std::vector<clique_t>
triangulate_graph(const mutk::RelationshipGraph &graph) {
    using LocalGraph = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS>;

    // Build up new graph
    LocalGraph local_graph(num_vertices(graph));
    for(auto v : make_vertex_range(graph)) {
        // add child vertex and link up all family members
        auto inv_range = make_inv_vertex_range(v, graph);
        for(auto w : inv_range) {
            add_edge(w, v, local_graph);
        }
        for(auto it = std::begin(inv_range); it != std::end(inv_range); ++it) {
            for(auto jt = std::next(it); jt != std::end(inv_range); ++jt) {
                add_edge(*it, *jt, local_graph);
            }
        }
    }

    auto fill_in_count = [](LocalGraph::vertex_descriptor v, const LocalGraph &g) {
        int fill = 0;
        auto adj_range = make_adj_vertex_range(v, g);
        for(auto it = std::begin(adj_range); it != std::end(adj_range); ++it) {
            for(auto jt = std::next(it); jt != std::end(adj_range); ++jt) {
                auto val = edge(*it, *jt, g);
                fill += !val.second;
            }
        }
        // return negative fill in because we are using a max heap
        return -fill;
    };
    using heap_value_t = std::pair<int, LocalGraph::vertex_descriptor>;
    using heap_t = boost::heap::d_ary_heap<heap_value_t,
        boost::heap::arity<2>, boost::heap::mutable_<true>>;

    heap_t priority_queue;

    std::vector<heap_t::handle_type> handles(num_vertices(graph));

    for(auto v : make_vertex_range(local_graph)) {
        int cost = fill_in_count(v, local_graph);
        auto handle = priority_queue.push({cost, v});
        handles[v] = handle;
    }

    std::vector<clique_t> elim_order;

    while(!priority_queue.empty()) {
        // chose next vertex to eliminate and remove it from the queue
        auto value = priority_queue.top();
        priority_queue.pop();

        // record the vertex
        auto v = value.second;
        auto & clique = elim_order.emplace_back();
        clique.push_back(v);
        // record neighbors
        for(auto n : make_adj_vertex_range(v, local_graph)) {
            clique.push_back(n);
        }
        // clear vertex
        clear_vertex(v, local_graph);

        // link up neighbors
        std::set<LocalGraph::vertex_descriptor> dirty_vertices;
        for(auto it = std::next(clique.begin()); it != clique.end(); ++it) {
            for(auto jt = std::next(it); jt != clique.end(); ++jt) {
                add_edge(*it, *jt, local_graph);
            }
            // keep track of vertices that need to be recalculated
            dirty_vertices.insert(*it);
            for(auto n : make_adj_vertex_range(*it, local_graph)) {
                dirty_vertices.insert(n);
            }
        }
        // Update the priority queue
        for(auto v : dirty_vertices) {
            (*handles[v]).first = fill_in_count(v, local_graph);
            priority_queue.update(handles[v]);
        }
    }

    return elim_order;
}

// LCOV_EXCL_START
TEST_CASE("triangulate_graph() identifies cliques") {
    using mutk::RelationshipGraph;
    using mutk::sample_id_t;

    {   
        RelationshipGraph graph(9);
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

        auto cliques = triangulate_graph(graph);

        REQUIRE(cliques.size() == 9);
        CHECK(cliques[0] == clique_t({8,6,7}));
        CHECK(cliques[1] == clique_t({5,4,7}));
        CHECK(cliques[2] == clique_t({2,3,6}));
        CHECK(cliques[3] == clique_t({7,4,6}));
        CHECK(cliques[4] == clique_t({6,3,4}));
        CHECK(cliques[5] == clique_t({4,0,1,3}));
        CHECK(cliques[6] == clique_t({3,0,1}));
        CHECK(cliques[7] == clique_t({1,0}));
        CHECK(cliques[8] == clique_t({0}));
    }
    {
        RelationshipGraph graph(7);
        add_edge(0,2,graph);
        add_edge(1,2,graph);
        add_edge(2,3,graph);
        add_edge(2,4,graph);
        add_edge(5,6,graph);

        auto cliques = triangulate_graph(graph);
        REQUIRE(cliques.size() == 7);
        CHECK(cliques[0] == clique_t({6,5}));
        CHECK(cliques[1] == clique_t({5}));
        CHECK(cliques[2] == clique_t({4,2}));
        CHECK(cliques[3] == clique_t({3,2}));
        CHECK(cliques[4] == clique_t({2,0,1}));
        CHECK(cliques[5] == clique_t({1,0}));
        CHECK(cliques[6] == clique_t({0}));
    }
}
// LCOV_EXCL_STOP

std::vector<component_t>
calculate_components(const mutk::RelationshipGraph &graph) {
    using mutk::variable_t;

    std::vector<component_t> components;

    for(auto v : make_vertex_range(graph)) {
        {
            // Add a potential for this single node
            auto & pot = components.emplace_back();
            pot.variables.push_back(variable_t(v));
            pot.edge_lengths.push_back(0.0f);
        }
        if(in_degree(v,graph) > 0) {
            // Add a potential for this family
            auto & pot = components.emplace_back();
            pot.variables.push_back(variable_t(v));
            pot.edge_lengths.push_back(0.0f);
            for(auto e : boost::make_iterator_range(in_edges(v,graph))) {
                pot.variables.push_back(variable_t(source(e, graph)));
                pot.edge_lengths.push_back(get(boost::edge_length, graph, e));
            }
        }
    }

    return components;
}