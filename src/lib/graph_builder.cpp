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

#include <boost/graph/copy.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/none_of.hpp>
#include <boost/algorithm/cxx11/one_of.hpp>
#include <boost/heap/d_ary_heap.hpp>

#include "mutk/graph_builder.hpp"

namespace {
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
}

using mutk::member_id_t;

using mutk::potential_t;

static mutk::relationship_graph::Graph
simplify_graph(mutk::relationship_graph::Graph &graph);

static std::vector<clique_t>
triangulate_graph(const mutk::relationship_graph::Graph &graph);

static std::vector<potential_t>
calculate_potentials(const mutk::relationship_graph::Graph &graph);

static mutk::relationship_graph::JunctionTree
create_junction_tree(const mutk::relationship_graph::Graph &graph,
    const std::vector<potential_t> &potentials,
    const std::vector<clique_t> &elimination_order);

// Convert `name` into a member id. If `name` is already registered, it
// return the registered id number. Otherwise add `name` to the registry.
member_id_t mutk::GraphBuilder::LookupName(const std::string &name) {
    auto it = map_member_name_to_id_.find(name);
    if(it != map_member_name_to_id_.end()) {
        return it->second;
    }

    member_id_t id{member_counter_++};
    map_member_name_to_id_.try_emplace(name, id);
    member_names_.push_back(name);

    // resize these vectors as well; add information later
    member_sexes_.push_back({});
    member_input_samples_.push_back({});
    return id;
}

member_id_t mutk::GraphBuilder::AddSingle(const std::string &name, const std::string &sex,
        const std::vector<std::string> &sample_names) {
    auto child_id = LookupName(name);
    member_sexes_[+child_id] = sex;
    member_input_samples_[+child_id] = sample_names;
    components_.push_back(component_t{{child_id}, {0.0f}});
    return child_id;
}

member_id_t mutk::GraphBuilder::AddPair(const std::string &name, const std::string &sex,
        const std::vector<std::string> &sample_names, const std::string &parent_name, float mutation_scale) {
    auto parent_id = LookupName(parent_name);
    auto child_id = LookupName(name);
    member_sexes_[+child_id] = sex;
    member_input_samples_[+child_id] = sample_names;
    components_.push_back(component_t{{child_id, parent_id}, {0.0f, mutation_scale}});
    return child_id;
}

member_id_t mutk::GraphBuilder::AddTrio(const std::string &name, const std::string &sex,
        const std::vector<std::string> &sample_names,
        const std::string &parent_name_a, float mutation_scale_a,
        const std::string &parent_name_b, float mutation_scale_b) {
    auto parent_id_a = LookupName(parent_name_a);
    auto parent_id_b = LookupName(parent_name_b);
    auto child_id = LookupName(name);
    member_sexes_[+child_id] = sex;
    member_input_samples_[+child_id] = sample_names;
    components_.push_back(component_t{{child_id, parent_id_a, parent_id_b},
        {0.0f, mutation_scale_a, mutation_scale_b}});
    return child_id;
}

void mutk::GraphBuilder::AddSamples(const std::vector<std::string> &sample_names) {
    sample_names_ = sample_names;
    map_sample_name_to_id_.clear();
    for(size_t i=0; i < sample_names.size(); ++i) {
        map_sample_name_to_id_.try_emplace(sample_names[i], sample_id_t{(int)i});
    }
}

void mutk::GraphBuilder::BuildGraph(const InheritanceModel &model, float mu) {
    using mutk::relationship_graph::Graph;
    static_assert(std::is_integral_v<Graph::vertex_descriptor>);

    auto initial_graph = CreateInitialGraph(model, mu);

    auto graph = simplify_graph(initial_graph);

    auto potentials = calculate_potentials(graph);

    auto cliques = triangulate_graph(graph);


}

mutk::relationship_graph::Graph mutk::GraphBuilder::CreateInitialGraph(const InheritanceModel &model, float mu) {
    using mutk::relationship_graph::Graph;

    Graph graph(member_counter_);

    auto data = get(boost::vertex_data, graph);
    auto ploidies = get(boost::vertex_ploidy, graph);
    auto labels = get(boost::vertex_label, graph);

    std::vector<InheritanceModel::chromosome_type_t> member_types;

    for(int i=0; i < member_counter_; ++i) {
        // Save node name
        labels[i] = member_names_[i];

        // Convert sex from strings to number
        auto it = model.map_name_to_type_.find(member_sexes_[i]);
        if(it == model.map_name_to_type_.end()) {
            throw std::invalid_argument("Member " + member_names_[i] +
                "has invalid sex: '" + member_sexes_[i] + "'.");
        }
        member_types.push_back(it->second);

        // Set ploidy
        ploidies[i] = mutk::relationship_graph::VertexPloidy{model.ploidies_[+it->second]};

        // Identify samples
        auto & member_data = member_data_samples_.emplace_back();
        for(auto & sample : member_input_samples_[i]) {
            auto it = map_sample_name_to_id_.find(sample);
            if(it != map_sample_name_to_id_.end()) {
                member_data.push_back(it->second);
            }
        }
        data[i] = member_data;
    }

    // Add edges
    for(const auto & component : components_) {
        // Lookup family structure in inheritance model
        std::vector<InheritanceModel::chromosome_type_t> pattern;
        for(const auto & member : component.members) {
            pattern.push_back(member_types[+member]);
        }
        auto pattern_it = std::find_if(model.patterns_.begin(), model.patterns_.end(), [&](const auto &val){
            return pattern == val.pattern;
        });
        if(pattern_it == model.patterns_.end()) {
            throw std::invalid_argument("Member " + member_names_[+component.members[0]] +
                "has invalid family pattern.");
        }

        // Connect up all members of this component
        auto child_id = component.members[0];
        for(size_t j = 1; j < component.members.size(); ++j) {
                member_id_t from = component.members[j];
                if(!pattern_it->discard[j]) {
                    add_edge(+from, +child_id, {mu*component.scales[j]}, graph);
                }
        }
    }

    return graph;
}

// Possible crosses
//  2 x 2 -> 2; aa x bb -> ab
//  2 x 2 -> 1; ?????????????
//  2 x 1 -> 1; aa x b. -> a.
//  2 x 1 -> 2; aa x b. -> ab
//  1 x 0 -> 1; a. x .. -> a.
//  1 x 0 -> 0; a. x .. -> ..
//  1 x 1 -> 2; a x b -> ab
//  1 x 1 -> 1; a x b -> a
//  2 -> 2; ab -> ab
//  1 -> 1; a -> a
//  0 -> 0; . -> .
//  2 -> 1; aa -> a

mutk::relationship_graph::Graph
simplify_graph(mutk::relationship_graph::Graph &graph) {
    using Graph = mutk::relationship_graph::Graph;
    using boost::make_iterator_range;

    // Topologically sort members in reverse order
    std::vector<Graph::vertex_descriptor> rev_topo_order;
    topological_sort(graph, std::back_inserter(rev_topo_order));
    auto topo_order = boost::adaptors::reverse(rev_topo_order);

    // Simplify the original graph
    auto lengths = get(boost::edge_length, graph);
    auto data = get(boost::vertex_data, graph);

    // Clear all leaf vertexes that do not have samples, starting from the tips
    for(auto v : rev_topo_order) {
        if(data[v].empty() && out_degree(v, graph) == 0) {
            clear_vertex(v, graph);
        }
    }

    // Remove surplus founders.
    for(auto && v : topo_order) {
        auto inv_range = boost::make_iterator_range(inv_adjacent_vertices(v, graph));
        // A vertex has no sibs if its parents only have one child
        // Or it has no parents.
        bool no_sibs = boost::algorithm::all_of(inv_range, [&](const auto & u) {
            return (data[u].empty() && degree(u,graph) == 1);
        });
        if(no_sibs) {
            // remove in edges
            clear_in_edges(v, graph);
        }
    }

    // try to bypass nodes with no data that have one out_edge
    for(auto && v : topo_order) {
        if(!data[v].empty() || in_degree(v,graph) == 0 || out_degree(v,graph) != 1) {
            continue;
        }
        auto in_edge_range = boost::make_iterator_range(in_edges(v,graph));
        auto out_edge = *out_edges(v,graph).first;
        auto child = target(out_edge, graph);
        // If the total in-degree of child and v is > 3 then we can't
        // simplify this node because child would have more than 2
        // in edges.
        if(in_degree(child,graph)+in_degree(v,graph) > 3) {
            continue;
        }
        auto len = lengths[out_edge];
        for(auto &&e : in_edge_range) {
            auto grand = source(e, graph);
            add_edge(grand, child, {len+lengths[e]}, graph);
        }
        clear_vertex(v, graph);
    }

    // Construct new graph in topological order
    Graph new_graph;

    auto vertex_all_in = get(boost::vertex_all, graph);
    auto vertex_all_out = get(boost::vertex_all, new_graph);
    auto edge_all_in = get(boost::edge_all, graph);
    auto edge_all_out = get(boost::edge_all, new_graph);

    std::vector<Graph::vertex_descriptor> old_to_new_vertex(num_vertices(graph));
    for(auto o : topo_order) {
        if(degree(o, graph) == 0 && data[o].size() == 0) {
            continue;
        }
        auto v = add_vertex(new_graph);
        old_to_new_vertex[o] = v;
        put(vertex_all_out, v, get(vertex_all_in, o));
    }

    for(auto e : boost::make_iterator_range(edges(graph))) {
        auto a = old_to_new_vertex[source(e, graph)];
        auto b = old_to_new_vertex[target(e, graph)];
        auto f = add_edge(a, b, new_graph);
        put(edge_all_out, f.first, get(edge_all_in, e));
    }

    return new_graph;
}

// LCOV_EXCL_START
TEST_CASE("simplify_graph() simplifies relationship graphs") {
    using mutk::relationship_graph::Graph;
    using mutk::sample_id_t;

    auto get_length = [](Graph::vertex_descriptor a, Graph::vertex_descriptor b,
        const Graph &g) -> float {
        auto e = edge(a,b,g);
        if(e.second) {
            return get(boost::edge_length, g, e.first);
        }
        return NAN;
    };

    auto get_name = [](Graph::vertex_descriptor a, const Graph &g) -> std::string {
        return get(boost::vertex_label, g, a);
    };

    {   // STANDARD TRIO
        Graph graph(3);
        add_edge(0,2,0.5f,graph);
        add_edge(1,2,1.0f,graph);

        auto labels = get(boost::vertex_label, graph);
        labels[0] = "A";
        labels[1] = "B";
        labels[2] = "C";

        auto data = get(boost::vertex_data, graph);
        data[0].push_back(sample_id_t{0});
        data[1].push_back(sample_id_t{1});
        data[2].push_back(sample_id_t{2});

        auto out_graph = simplify_graph(graph);

        CHECK(num_vertices(out_graph) == 3);
        CHECK(get_name(0,out_graph) == "B");
        CHECK(get_name(1,out_graph) == "A");
        CHECK(get_name(2,out_graph) == "C");

        CHECK(num_edges(out_graph) == 2);
        CHECK(get_length(0,2,out_graph) == 1.0f);
        CHECK(get_length(1,2,out_graph) == 0.5f);
    }
    {   // TRIO WITH NO CHILD DATA
        Graph graph(3);
        add_edge(0,2,0.5f,graph);
        add_edge(1,2,1.0f,graph);

        auto labels = get(boost::vertex_label, graph);
        labels[0] = "A";
        labels[1] = "B";
        labels[2] = "C";

        auto data = get(boost::vertex_data, graph);
        data[0].push_back(sample_id_t{0});
        data[1].push_back(sample_id_t{1});

        auto out_graph = simplify_graph(graph);
        CHECK(num_vertices(out_graph) == 2);
        CHECK(get_name(0,out_graph) == "B");
        CHECK(get_name(1,out_graph) == "A");

        CHECK(num_edges(out_graph) == 0);
    }
    {   // TRIO WITH NO PARENT DATA
        Graph graph(3);
        add_edge(0,2,0.5f,graph);
        add_edge(1,2,1.0f,graph);

        auto labels = get(boost::vertex_label, graph);
        labels[0] = "A";
        labels[1] = "B";
        labels[2] = "C";

        auto data = get(boost::vertex_data, graph);
        data[2].push_back(sample_id_t{2});

        auto out_graph = simplify_graph(graph);

        CHECK(num_vertices(out_graph) == 1);
        CHECK(get_name(0,out_graph) == "C");

        CHECK(num_edges(out_graph) == 0);
    }
    {   // TRIO WITH EXTRA NODE
        Graph graph(4);
        add_edge(0,2,0.5f,graph);
        add_edge(1,2,1.0f,graph);
        add_edge(2,3,0.1f,graph);

        auto labels = get(boost::vertex_label, graph);
        labels[0] = "A";
        labels[1] = "B";
        labels[2] = "X";
        labels[3] = "C";

        auto data = get(boost::vertex_data, graph);
        data[0].push_back(sample_id_t{0});
        data[1].push_back(sample_id_t{1});
        data[3].push_back(sample_id_t{2});

        auto out_graph = simplify_graph(graph);

        CHECK(num_vertices(out_graph) == 3);
        CHECK(get_name(0,out_graph) == "B");
        CHECK(get_name(1,out_graph) == "A");
        CHECK(get_name(2,out_graph) == "C");

        CHECK(num_edges(out_graph) == 2);
        CHECK(get_length(0,2,out_graph) == 1.1f);
        CHECK(get_length(1,2,out_graph) == 0.6f);
    }
    {   // MONOZYGOTIC TWINS
        Graph graph(5);
        add_edge(0,2,0.5f,graph);
        add_edge(1,2,1.0f,graph);
        add_edge(2,3,0.1f,graph);
        add_edge(2,4,0.2f,graph);

        auto labels = get(boost::vertex_label, graph);
        labels[0] = "A";
        labels[1] = "B";
        labels[2] = "X";
        labels[3] = "C";
        labels[4] = "D";

        auto data = get(boost::vertex_data, graph);
        data[0].push_back(sample_id_t{0});
        data[1].push_back(sample_id_t{1});
        data[3].push_back(sample_id_t{2});
        data[4].push_back(sample_id_t{3});

        auto out_graph = simplify_graph(graph);

        CHECK(num_vertices(out_graph) == 5);
        CHECK(get_name(0,out_graph) == "B");
        CHECK(get_name(1,out_graph) == "A");
        CHECK(get_name(2,out_graph) == "X");
        CHECK(get_name(3,out_graph) == "D");
        CHECK(get_name(4,out_graph) == "C");

        CHECK(num_edges(out_graph) == 4);
        CHECK(get_length(0,2,out_graph) == 1.0f);
        CHECK(get_length(1,2,out_graph) == 0.5f);
        CHECK(get_length(2,3,out_graph) == 0.2f);
        CHECK(get_length(2,4,out_graph) == 0.1f);
    }
}
// LCOV_EXCL_STOP


// Triangulate a graph where the vertices are in topological order
//
// Almond and Kong (1991) Optimality Issues in Constructing a Markov Tree from Graphical Models.
//     Research Report 329. University of Chicago, Dept. of Statistics
static std::vector<clique_t>
triangulate_graph(const mutk::relationship_graph::Graph &graph) {
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
    using mutk::relationship_graph::Graph;
    using mutk::sample_id_t;

    {   
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
        Graph graph(7);
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

std::vector<potential_t>
calculate_potentials(const mutk::relationship_graph::Graph &graph) {
    std::vector<potential_t> potentials;

    for(auto v : make_vertex_range(graph)) {
        {
            // Add a potential for this single node
            auto & pot = potentials.emplace_back();
            pot.variables.emplace_back(v);
            pot.edge_lengths.push_back(0.0f);
        }
        if(in_degree(v,graph) > 0) {
            // Add a potential for this family
            auto & pot = potentials.emplace_back();
            pot.variables.emplace_back(v);
            pot.edge_lengths.push_back(0.0f);
            for(auto e : boost::make_iterator_range(in_edges(v,graph))) {
                pot.variables.emplace_back(source(e, graph));
                pot.edge_lengths.push_back(get(boost::edge_length, graph, e));
            }
        }
    }

    return potentials;
}
