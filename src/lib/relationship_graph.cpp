/*
# Copyright (c) 2014-2020 Reed A. Cartwright <reed@cartwright.ht>
# Copyright (c) 2016 Steven H. Wu
# Copyright (c) 2016 Kael Dai
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

#include <mutk/relationship_graph.hpp>
#include <mutk/detail/graph.hpp>

#include <cmath>
#include <queue>
#include <algorithm>
#include <unordered_set>

#include <boost/range/algorithm/find.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/heap/d_ary_heap.hpp>
#include <boost/range/algorithm/min_element.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/range/algorithm/fill.hpp>

namespace {

using Sex = mutk::Pedigree::Sex;
using IndyId = mutk::IndyId;

using samples_t = mutk::RelationshipGraph::samples_t;

using InheritanceModel = mutk::InheritanceModel;
using Pedigree = mutk::Pedigree;

namespace pedigree_graph = mutk::detail::pedigree_graph;
namespace junction_tree = mutk::detail::junction_tree;
using neighbors_t = junction_tree::neighbors_t;

using VertexType = mutk::detail::VertexType;
using EdgeType = mutk::detail::EdgeType;

void construct_pedigree_graph(pedigree_graph::Graph &graph,
    const Pedigree &pedigree, const samples_t& known_samples,
    bool normalize_somatic_trees);

void update_edge_lengths(pedigree_graph::Graph &graph, double mu_meiotic, double mu_somatic);
void simplify(pedigree_graph::Graph &graph);
void prune(pedigree_graph::Graph &graph, InheritanceModel model);

pedigree_graph::Graph finalize(const pedigree_graph::Graph &input);

std::vector<mutk::potential_t>
create_potentials(const pedigree_graph::Graph &graph);

std::pair<std::vector<pedigree_graph::vertex_t>,
std::vector<neighbors_t>>
triangulate_graph(const pedigree_graph::Graph &graph);

junction_tree::Graph create_junction_tree(
    const std::vector<pedigree_graph::vertex_t> &elim_order,
    const std::vector<neighbors_t> &neighbors);

} // namespace

const std::map<std::string, InheritanceModel> mutk::detail::CHR_MODEL_MAP {
    {"autosomal", InheritanceModel::Autosomal},
    {"maternal", InheritanceModel::Maternal},
    {"paternal", InheritanceModel::Paternal},
    {"x-linked", InheritanceModel::XLinked},
    {"y-linked", InheritanceModel::YLinked},
    {"w-linked", InheritanceModel::WLinked},
    {"z-linked", InheritanceModel::ZLinked},
    {"mitochondrial", InheritanceModel::Maternal},
    {"xlinked", InheritanceModel::XLinked},
    {"ylinked", InheritanceModel::YLinked},
    {"wlinked", InheritanceModel::WLinked},
    {"zlinked", InheritanceModel::ZLinked},
};

float mutk::RelationshipGraph::PeelForward(workspace_t &work) const {
    assert(work.stack.size() >= stack_size_);

    for(auto &&peeler : peelers_) {
        peeler->Forward(work);
    }
    float log_ret = work.scale;

    for(auto r : roots_) {
        float s = xt::sum(work.stack[r])();
        log_ret += std::log(s);
    }
    return log_ret;
}

void mutk::RelationshipGraph::ConstructGraph(const Pedigree& pedigree,
        const samples_t& known_samples, InheritanceModel model,
        double mu, double mu_somatic,
        bool normalize_somatic_trees) {
    using namespace std;

    inheritance_model_ = model;

    // Construct a boost::graph of the pedigree and somatic information
    pedigree_graph::Graph graph;

    construct_pedigree_graph(graph, pedigree, known_samples, normalize_somatic_trees);
 
    // Multiply edge lengths by mutation rates
    update_edge_lengths(graph, mu, mu_somatic);

    // Remove edges that are non-informative
    simplify(graph);

    // Prune pedigree
    prune(graph, inheritance_model_);

    // Sort and eliminate cleared vertices
    graph_ = finalize(graph);

    // Identify founders, non-founders, and samples
    founders_.first = 0;
    for(founders_.second=0; founders_.second < num_vertices(graph_); ++founders_.second) {
        if(in_degree(founders_.second,graph_) > 0) {
            break;
        }
    }
    descendants_.first = founders_.second;
    for(descendants_.second=0; descendants_.second < num_vertices(graph_); ++descendants_.second) {
        if(out_degree(descendants_.second,graph_) == 0) {
            break;
        }
    }
    samples_.first = descendants_.second;
    samples_.second = num_vertices(graph_);
}

// LCOV_EXCL_START
namespace {
class Test_ConstructGraph : public mutk::RelationshipGraph {
TEST_CASE_CLASS("RelationshipGraph-ConstructGraph") {
    SUBCASE("InheritanceModel::Autosomal") {
        const char ped[] = 
            "##PEDNG v0.1\n"
            "A\t.\t.\tM\t=\n"
            "B\t.\t.\tF\t=\n"
            "C\tA\tB\tM\t=\n";

        std::vector<const char*> known_samples = {"A","B","C"};

        Pedigree pedigree;
        REQUIRE_NOTHROW(pedigree = Pedigree::parse_text(ped));

        Test_ConstructGraph relationship_graph;

        REQUIRE_NOTHROW(relationship_graph.ConstructGraph(pedigree,
            known_samples, mutk::InheritanceModel::Autosomal,
            1e-6, 1e-6, false));

        auto graph = relationship_graph.graph();

                auto labels = get(boost::vertex_label, graph);
        auto ploidies = get(boost::vertex_ploidy, graph);
        auto types = get(boost::vertex_type, graph);

        CHECK(labels[0] == "B/z");
        CHECK(labels[1] == "A/z");
        CHECK(labels[2] == "B");
        CHECK(labels[3] == "C");
        CHECK(labels[4] == "A");

        CHECK(ploidies[0] == 2);
        CHECK(ploidies[1] == 2);
        CHECK(ploidies[2] == 2);
        CHECK(ploidies[3] == 2);
        CHECK(ploidies[4] == 2);

        CHECK(types[0] == VertexType::Germline);
        CHECK(types[1] == VertexType::Germline);
        CHECK(types[2] == VertexType::Sample);
        CHECK(types[3] == VertexType::Sample);
        CHECK(types[4] == VertexType::Sample);

        auto adj0 = boost::make_iterator_range(adjacent_vertices(0, graph));
        auto inv0 = boost::make_iterator_range(inv_adjacent_vertices(0, graph));
        REQUIRE(adj0.size() == 2);
        REQUIRE(inv0.size() == 0);
        CHECK(adj0[0] == 2);
        CHECK(adj0[1] == 3);

        auto adj1 = boost::make_iterator_range(adjacent_vertices(1, graph));
        auto inv1 = boost::make_iterator_range(inv_adjacent_vertices(1, graph));
        REQUIRE(adj1.size() == 2);
        REQUIRE(inv1.size() == 0);
        CHECK(adj1[0] == 4);
        CHECK(adj1[1] == 3);

        auto adj2 = boost::make_iterator_range(adjacent_vertices(2, graph));
        auto inv2 = boost::make_iterator_range(inv_adjacent_vertices(2, graph));
        REQUIRE(adj2.size() == 0);
        REQUIRE(inv2.size() == 1);
        CHECK(inv2[0] == 0);

        auto adj3 = boost::make_iterator_range(adjacent_vertices(3, graph));
        auto inv3 = boost::make_iterator_range(inv_adjacent_vertices(3, graph));
        REQUIRE(adj3.size() == 0);
        REQUIRE(inv3.size() == 2);
        CHECK(inv3[0] == 1);
        CHECK(inv3[1] == 0);

        auto adj4 = boost::make_iterator_range(adjacent_vertices(4, graph));
        auto inv4 = boost::make_iterator_range(inv_adjacent_vertices(4, graph));
        REQUIRE(adj4.size() == 0);
        REQUIRE(inv4.size() == 1);
        CHECK(inv3[0] == 1);

        CHECK(relationship_graph.founders_.first == 0);
        CHECK(relationship_graph.founders_.second == 2);
        CHECK(relationship_graph.descendants_.first == 2);
        CHECK(relationship_graph.descendants_.second == 2);
        CHECK(relationship_graph.samples_.first == 2);
        CHECK(relationship_graph.samples_.second == 5);
    }
}
};
} //namespace anon
// LCOV_EXCL_STOP

void mutk::RelationshipGraph::ConstructPeeler() {
    potentials_ = create_potentials(graph_);

    auto elim_steps = triangulate_graph(graph_);

    junction_tree_ = create_junction_tree(elim_steps.first, elim_steps.second);

    // DEBUG: Print Junction Tree
    // std::cerr << std::endl;    
    // for(auto &&e : boost::make_iterator_range(edges(junction_tree_))) {
    //     auto v = source(e,junction_tree_);
    //     bool b = false;
    //     for(auto w : get(boost::vertex_label, junction_tree_, v)) {
    //         if(b) {
    //             std::cerr << ", ";
    //         }
    //         b = true;
    //         std::cerr << get(boost::vertex_label, graph_, w);
    //     }

    //     std::cerr << " --- ";
    //     v = target(e,junction_tree_);
    //     b = false;
    //     for(auto w : get(boost::vertex_label, junction_tree_, v)) {
    //         if(b) {
    //             std::cerr << ", ";
    //         }
    //         b = true;
    //         std::cerr << get(boost::vertex_label, graph_, w);
    //     }
    //     std::cerr << std::endl;
    // }

    // build a container that aligns a set of graph vertices with a node in the
    // junction tree
    boost::container::flat_map<junction_tree::neighbors_t,junction_tree::vertex_t>
        map_labels_to_jctnodes;
    auto jctnode_labels = get(boost::vertex_label, junction_tree_);
    for(junction_tree::vertex_t v = 0; v < num_vertices(junction_tree_); ++v) {
        map_labels_to_jctnodes[jctnode_labels[v]] = v;
    }
    // construct the set of graph vertices for each potential, and lookup the
    // junction tree node corresponding to the set.
    std::vector<size_t> jctnode_to_potential(num_vertices(junction_tree_), -1);
    for(size_t i = 0; i < potentials_.size(); ++i) {
        neighbors_t label;
        for(auto &&a : potentials_[i].parents) {
            label.insert(+a.first);
        }
        label.insert(+potentials_[i].child);
        junction_tree::vertex_t v = map_labels_to_jctnodes[label];
        // check for duplicated potentials
        assert(jctnode_to_potential[v] == -1);
        jctnode_to_potential[v] = i;
    }
    // For nodes missing potentials add some unit ones.
    for(junction_tree::vertex_t v = 0; v < jctnode_to_potential.size(); ++v) {
        if(jctnode_to_potential[v] != -1) {
            continue;
        }
        potential_t unit;
        unit.child = IndyId(-1);
        unit.type = PotentialType::Unit;
        auto label = jctnode_labels[v];
        for(auto &&a : label) {
            unit.parents.emplace_back(IndyId(a), 0.0f);
        }
        potentials_.push_back(unit);
        jctnode_to_potential[v] = potentials_.size()-1;
    }

    // sort the labels of each jctnode according to elimination rank
    // elimination_rank_[id] = at what step `id` is eliminated
    elimination_rank_.resize(num_vertices(graph_));
    for(size_t n = 0; n < elim_steps.first.size(); ++n) {
        elimination_rank_[elim_steps.first[n]] = n;
    }
    
    auto ploidies = get(boost::vertex_ploidy, graph_);
    using tensor_label_t = std::vector<PeelingVertex::label_t>;
    std::vector<tensor_label_t> tensor_labels(num_vertices(junction_tree_));
    auto jctnode_range = boost::make_iterator_range(vertices(junction_tree_));

    // sort the vertices in jctnode by elimination rank
    // vertices that are eliminated later go to the front.
    // use a lambda function for code reuse
    auto elim_cmp = [](auto a, auto b) -> bool {
        return a > b;
    };

    for(auto &&v : jctnode_range) {
        // sort the vertices at this jctnode by elimination rank
        const neighbors_t &label = jctnode_labels[v];
        tensor_label_t tensor_lab(label.begin(), label.end());
        boost::range::sort(tensor_lab, [&](auto a, auto b) {
            return elim_cmp(elimination_rank_[a], elimination_rank_[b]);
        });
        tensor_labels[v] = tensor_lab;

        auto &pot = potentials_[jctnode_to_potential[v]];

        if(pot.type == PotentialType::Unit) {
            // for unit potentials, use axes to store the node shape as ploidies
            pot.axes.clear();
            for(auto && lab : tensor_lab) {
                pot.axes.push_back(ploidies[lab]);
            }
            continue;
        }

        // Transition matrices when created have axes child,parent,parent
        // identify the proper shuffle axes.
        std::vector<std::pair<std::size_t,std::size_t>> labeled_axes;
        labeled_axes.emplace_back(elimination_rank_[+pot.child],0);
        for(size_t i = 0; i < pot.parents.size(); ++i) {
            labeled_axes.emplace_back(+pot.parents[i].first,i+1);
        }
        // sort the labeled axes so that labels that are eliminated
        // later come earlier.
        boost::range::sort(labeled_axes, [&](auto a, auto b) {
            return elim_cmp(a.first, b.first);
        });

        pot.axes.resize(labeled_axes.size());
        for(size_t i = 0; i < labeled_axes.size(); ++i) {
            pot.axes[i] = labeled_axes[i].second;
        }
    }

    // The first part of the stack refers to potentials
    // Counter will track all the additional messages
    // that we need to store during peeling
    size_t counter = potentials_.size();

    // Record the output messages of every node
    std::vector<size_t> output_index(num_vertices(junction_tree_), -1);

    peeling_ops_.clear();

    // iterate the junction tree in reverse order
    auto rev_jctnode_range = boost::adaptors::reverse(jctnode_range);
    for(auto &&v : rev_jctnode_range) {
        peeling_op_t op;
        // Setup local data
        op.local.register_id = jctnode_to_potential[v];
        op.local.axes.assign(tensor_labels[v].begin(), tensor_labels[v].end());

        // Setup inputs and output connections
        auto adj_range = boost::make_iterator_range(adjacent_vertices(v, junction_tree_));
        for(auto &&w : adj_range) {
            assert(w != v);
            if(w < v) {
                assert(output_index[v] == -1);
                op.output.register_id = counter;
                op.output.axes.assign(tensor_labels[w].begin(), tensor_labels[w].end());
                output_index[v] = counter++;
            } else {
                peeling_op_t::data_t dat;
                dat.register_id = output_index[w];
                dat.axes.assign(tensor_labels[w].begin(), tensor_labels[w].end());
                op.inputs.push_back(dat);
            }
        }
        // If no output has been found, we have a root peeler
        if(output_index[v] == -1) {
            op.output.register_id = counter;
            op.output.axes.clear();
            output_index[v] = counter++;
        }
        peeling_ops_.push_back(op);
    }
    stack_size_ = counter;

    // Construct peelers
    peelers_.clear();
    for(const auto &op : peeling_ops_) {
        peeling_t peeler = std::make_unique<GeneralPeelingVertex>(op.local.axes);
        peeler->AddLocal(op.local.axes, op.local.register_id);
        for(const auto & i : op.inputs) {
            peeler->AddInput(i.axes, i.register_id);
        }
        peeler->AddOutput(op.output.axes, op.output.register_id);
        if(op.output.axes.empty()) {
            roots_.push_back(op.output.register_id);            
        }
        peelers_.push_back(std::move(peeler));
    }
}

// LCOV_EXCL_START
namespace {
class Test_ConstructPeeler : public mutk::RelationshipGraph {
TEST_CASE_CLASS("RelationshipGraph-ConstructPeeler Single") {
    using PotentialType = mutk::PotentialType;
    const char ped[] = 
        "##PEDNG v0.1\n"
        "A\t.\t.\tF\t=\n";
    std::vector<const char*> known_samples = {"A"};

    Pedigree pedigree;
    REQUIRE_NOTHROW(pedigree = Pedigree::parse_text(ped));

    SUBCASE("InheritanceModel::Autosomal") {
        Test_ConstructPeeler relationship_graph;
        REQUIRE_NOTHROW(relationship_graph.ConstructGraph(pedigree,
            known_samples, mutk::InheritanceModel::Autosomal,
            1e-6, 1e-6, false));
        REQUIRE_NOTHROW(relationship_graph.ConstructPeeler());

        // Checking Potentials
        auto &potentials = relationship_graph.potentials_;
        REQUIRE(potentials.size() == 1+1+1);
        CHECK(potentials[0].type == PotentialType::LikelihoodDiploid);
        CHECK(+potentials[0].child == 1);
        CHECK(potentials[0].parents.size() == 0);

        CHECK(potentials[1].type == PotentialType::FounderDiploid);
        CHECK(+potentials[1].child == 0);
        CHECK(potentials[1].parents.size() == 0);

        CHECK(potentials[2].type == PotentialType::CloneDiploid);
        CHECK(+potentials[2].child == 1);
        REQUIRE(potentials[2].parents.size() == 1);
        CHECK(+potentials[2].parents[0].first == 0);
        CHECK(potentials[2].parents[0].second == 1e-6f);

        // Elimination Rank
        auto &elimination_rank = relationship_graph.elimination_rank_;
        REQUIRE(elimination_rank.size() == 2);
        CHECK(elimination_rank[0] == 1);
        CHECK(elimination_rank[1] == 0);

        // Checking Peelers
        auto &peelers = relationship_graph.peelers_;
        REQUIRE(peelers.size() == 3);
        mutk::GeneralPeelingVertex *peel0 = dynamic_cast<mutk::GeneralPeelingVertex*>(peelers[0].get());
        CHECK(peel0->local_data().index == 0);
        CHECK(peel0->output_data().index == 3);
        CHECK(peel0->input_data().empty());
        mutk::GeneralPeelingVertex *peel1 = dynamic_cast<mutk::GeneralPeelingVertex*>(peelers[1].get());
        CHECK(peel1->local_data().index == 2);
        CHECK(peel1->output_data().index == 4);
        REQUIRE(peel1->input_data().size() == 1);
        CHECK(peel1->input_data()[0].index == 3);
        mutk::GeneralPeelingVertex *peel2 = dynamic_cast<mutk::GeneralPeelingVertex*>(peelers[2].get());
        CHECK(peel2->local_data().index == 1);
        CHECK(peel2->output_data().index == 5);
        REQUIRE(peel2->input_data().size() == 1);
        CHECK(peel2->input_data()[0].index == 4);

        // Checking Root
        REQUIRE(relationship_graph.roots_.size() == 1);
        CHECK(relationship_graph.roots_[0] == 5);

        // Checking Product
        auto work = relationship_graph.CreateWorkspace();
        work.scale = 10.0;
        work.stack[0] = {1, 0.1, 0.001};
        work.stack[1] = {1-1e-4-1e-8, 1e-4, 1e-8};
        work.stack[2] = {{0.998001, 0.000999, 1e-06},
                         {0.001998, 0.998002, 0.001998},
                         {1e-06, 0.000999, 0.998001}};

        CHECK(relationship_graph.PeelForward(work) == doctest::Approx(10+std::log(0.9981111)));
    }
}

TEST_CASE_CLASS("RelationshipGraph-ConstructPeeler Trio") {
    using PotentialType = mutk::PotentialType;

    const char ped[] = 
        "##PEDNG v0.1\n"
        "A\t.\t.\tM\t=\n"
        "B\t.\t.\tF\t=\n"
        "C\tA\tB\tM\t=\n";
    std::vector<const char*> known_samples = {"A","B","C"};

    Pedigree pedigree;
    REQUIRE_NOTHROW(pedigree = Pedigree::parse_text(ped));

    SUBCASE("InheritanceModel::Autosomal") {
        Test_ConstructPeeler relationship_graph;
        REQUIRE_NOTHROW(relationship_graph.ConstructGraph(pedigree,
            known_samples, mutk::InheritanceModel::Autosomal,
            1e-6, 1e-6, false));
        REQUIRE_NOTHROW(relationship_graph.ConstructPeeler());

        SUBCASE("Verify Potentials") {
            using v = std::vector<int>;
            auto &potentials = relationship_graph.potentials_;

            REQUIRE(potentials.size() == 3+2+3+1);

            CHECK(potentials[0].type == PotentialType::LikelihoodDiploid);
            CHECK(+potentials[0].child == 2);
            CHECK(potentials[0].parents.size() == 0);

            CHECK(potentials[1].type == PotentialType::LikelihoodDiploid);
            CHECK(+potentials[1].child == 3);
            CHECK(potentials[1].parents.size() == 0);

            CHECK(potentials[2].type == PotentialType::LikelihoodDiploid);
            CHECK(+potentials[2].child == 4);
            CHECK(potentials[2].parents.size() == 0);

            CHECK(potentials[3].type == PotentialType::FounderDiploid);
            CHECK(+potentials[3].child == 0);
            CHECK(potentials[3].parents.size() == 0);

            CHECK(potentials[4].type == PotentialType::FounderDiploid);
            CHECK(+potentials[4].child == 1);
            CHECK(potentials[4].parents.size() == 0);

            CHECK(potentials[5].type == PotentialType::CloneDiploid);
            CHECK(+potentials[5].child == 2);
            REQUIRE(potentials[5].parents.size() == 1);
            CHECK(+potentials[5].parents[0].first == 0);
            CHECK(potentials[5].parents[0].second == 1e-6f);

            CHECK(potentials[6].type == PotentialType::ChildDiploidDiploid);
            CHECK(+potentials[6].child == 3);
            REQUIRE(potentials[6].parents.size() == 2);
            CHECK(+potentials[6].parents[0].first == 1);
            CHECK(potentials[6].parents[0].second == 2e-6f);
            CHECK(+potentials[6].parents[1].first == 0);
            CHECK(potentials[6].parents[1].second == 2e-6f);
            CHECK_EQ_RANGES(potentials[6].axes, v({0,1,2}));

            CHECK(potentials[7].type == PotentialType::CloneDiploid);
            CHECK(+potentials[7].child == 4);
            REQUIRE(potentials[7].parents.size() == 1);
            CHECK(+potentials[7].parents[0].first == 1);
            CHECK(potentials[7].parents[0].second == 1e-6f);

            CHECK(potentials[8].type == PotentialType::Unit);
            CHECK(+potentials[8].child == -1);
            REQUIRE(potentials[8].parents.size() == 2);
            CHECK(+potentials[8].parents[0].first == 0);
            CHECK(potentials[8].parents[0].second == 0.0f);
            CHECK(+potentials[8].parents[1].first == 1);
            CHECK(potentials[8].parents[1].second == 0.0f);
        }
        SUBCASE("Verify Peelers") {
            using v = std::vector<int>;
            CHECK_EQ_RANGES(relationship_graph.elimination_rank_, v({4,3,2,1,0}));
            CHECK_EQ_RANGES(relationship_graph.roots_, v({17}));
            
            auto &peelers = relationship_graph.peeling_ops_;

            REQUIRE(peelers.size() == 9);
            {
                std::vector<size_t> observed(9);
                std::transform(peelers.begin(), peelers.end(), observed.begin(), [&](auto & x){
                    return x.local.register_id;
                });
                CHECK_EQ_RANGES(observed, v({2,7,1,6,0,5,4,8,3}));
            } {
                std::vector<size_t> observed(9);
                std::transform(peelers.begin(), peelers.end(), observed.begin(), [&](auto & x){
                    return x.output.register_id;
                });
                CHECK_EQ_RANGES(observed, v({9,10,11,12,13,14,15,16,17}));
            } {
                auto input_ids = [&](auto & a) {
                    std::vector<size_t> observed;
                    std::transform(a.inputs.begin(), a.inputs.end(),
                    std::back_inserter(observed), [&](auto & x){
                        return x.register_id;
                    });
                    return observed;
                };
                CHECK_EQ_RANGES(input_ids(peelers[0]), v({}));
                CHECK_EQ_RANGES(input_ids(peelers[1]), v({9}));
                CHECK_EQ_RANGES(input_ids(peelers[2]), v({}));
                CHECK_EQ_RANGES(input_ids(peelers[3]), v({11}));
                CHECK_EQ_RANGES(input_ids(peelers[4]), v({}));
                CHECK_EQ_RANGES(input_ids(peelers[5]), v({13}));
                CHECK_EQ_RANGES(input_ids(peelers[6]), v({10}));
                CHECK_EQ_RANGES(input_ids(peelers[7]), v({15,12}));
                CHECK_EQ_RANGES(input_ids(peelers[8]), v({16,14}));
            }

        }
        SUBCASE("Verify Results") {
            //auto &potentials = relationship_graph.potentials_;
            auto work = relationship_graph.CreateWorkspace();
            work.scale = 0;
            // Likelihoods
            work.stack[0] = {1, 0, 0};    // Mom
            work.stack[1] = {0, 1, 1e-5}; // Child
            work.stack[2] = {1, 1e-4, 0}; // Dad
            // Founders
            work.stack[3] = {1-1e-4-1e-8, 1e-4, 1e-8};
            work.stack[4] = {1-1e-4-1e-8, 1e-4, 1e-8};
            // Transitions
            work.stack[5] = {{0.998001, 0.000999, 1e-06},
                             {0.001998, 0.998002, 0.001998},
                             {1e-06, 0.000999, 0.998001}};
            mutk::tensor_t temp = {{{0.99, 0.495, 0}, {0.01, 0.5, 0.99}, {0, 0.005, 0.01}},
                             {{0.5, 0.25, 0}, {0.5, 0.5, 0.5}, {0, 0.25, 0.5}},
                             {{0.01, 0.005, 0}, {0.99, 0.5, 0.01}, {0, 0.495, 0.99}}};
            work.stack[6] = xt::transpose(temp, {0,2,1});

            work.stack[7] = {{1, 0, 0},{0, 1, 0},{0, 0, 1}};
            // Unit potential
            work.stack[8] = xt::ones<float>({3,3});

            CHECK(relationship_graph.PeelForward(work) == doctest::Approx(std::log(0.009978069)));
        }
    }
}
};
} // anon namespace
// LCOV_EXCL_STOP

mutk::RelationshipGraph::workspace_t
mutk::RelationshipGraph::CreateWorkspace() const {
    workspace_t work;
    work.stack.resize(stack_size_);
    return work;    
}

mutk::RelationshipGraph::samples_t
mutk::RelationshipGraph::SampleNames() const {
    samples_t ret;

    auto vertex_range = boost::make_iterator_range(vertices(graph_));
    auto labels = get(boost::vertex_label, graph_);
    auto types = get(boost::vertex_type, graph_);

    for(auto v : vertex_range) {
        if(in_degree(v,graph_) > 0 && types[v] == VertexType::Sample) {
            ret.push_back(labels[v].c_str());
        }
    }
    return ret;
}

void mutk::RelationshipGraph::PrintGraph(std::ostream &os) const {
    auto vertex_range = boost::make_iterator_range(vertices(graph_));

    auto labels = get(boost::vertex_label, graph_);
    auto sexes = get(boost::vertex_sex, graph_);
    auto ploidies = get(boost::vertex_ploidy, graph_);
    auto types = get(boost::vertex_type, graph_);


    auto sex = [&](mutk::Pedigree::Sex s) {
        switch(s) {
        case Sex::Autosomal:
            return "autosomal";
        case Sex::Male:
            return "male";
        case Sex::Female:
            return "female";
        default:
            break;
        };
        return "unknown";
    };

    auto print_vertex = [&](auto &&v) {
        os << "  " << labels[v] << ":\n"
           << "    sex: " << sex(sexes[v]) << "\n"
           << "    ploidy: " << ploidies[v] << "\n";
        auto [ei,ee] = in_edges(v, graph_);
        if(ei == ee) {
            return;
        }
        os << "    origin:\n";
        for(auto it = ei; it != ee; ++it) {
            os << "      - label:  " << labels[source(*it,graph_)] << "\n"
               << "        length: " << get(boost::edge_length, graph_, *it) << "\n"
               << "        sex:    " << sex(sexes[source(*it,graph_)]) << "\n";
        }
    };

    os << "%YAML 1.2\n---\n";

    os << "founding:\n";
    for(auto v : vertex_range) {
        if(in_degree(v,graph_) == 0) {
            // add samples
            print_vertex(v);
        }
    }

    os << "\ngermline:\n";
    for(auto v : vertex_range) {
        if(in_degree(v,graph_) > 0 && types[v] == VertexType::Germline) {
            // add samples
            print_vertex(v);
        }
    }

    os << "\nsomatic:\n";
    for(auto v : vertex_range) {
        if(in_degree(v,graph_) > 0 && types[v] == VertexType::Somatic) {
            // add samples
            print_vertex(v);
        }
    }

    os << "\nsample:\n";
    for(auto v : vertex_range) {
        if(in_degree(v,graph_) > 0 && types[v] == VertexType::Sample) {
            // add samples
            print_vertex(v);
        }
    }
}

namespace {

void prune_autosomal(pedigree_graph::Graph &graph);
void prune_ylinked(pedigree_graph::Graph &graph);
void prune_xlinked(pedigree_graph::Graph &graph);
void prune_wlinked(pedigree_graph::Graph &graph);
void prune_zlinked(pedigree_graph::Graph &graph);
void prune_maternal(pedigree_graph::Graph &graph);
void prune_paternal(pedigree_graph::Graph &graph);

void prune(pedigree_graph::Graph &graph, InheritanceModel model) {
    switch(model) {
    case InheritanceModel::Autosomal:
        return prune_autosomal(graph);
    case InheritanceModel::YLinked:
        return prune_ylinked(graph);
    case InheritanceModel::XLinked:
        return prune_xlinked(graph);
    case InheritanceModel::WLinked:
        return prune_wlinked(graph);
    case InheritanceModel::ZLinked:
        return prune_zlinked(graph);
    case InheritanceModel::Maternal:
        return prune_maternal(graph);
    case InheritanceModel::Paternal:
        return prune_paternal(graph);
    default:
        throw std::invalid_argument("Selected inheritance model invalid or not implemented yet.");
    }
}

void prune_autosomal(pedigree_graph::Graph &graph) {
    /*do nothing*/;
}

void prune_ylinked(pedigree_graph::Graph &graph) {
    auto sexes  = get(boost::vertex_sex, graph);
    auto ploidies  = get(boost::vertex_ploidy, graph);

    auto is_x = [&](auto &&e) {
        if(!(get(boost::edge_type, graph, e) & mutk::detail::GERM_EDGE)) {
            return false;
        }
        auto a = source(e, graph);
        auto b = target(e, graph);
        return (sexes[a] == Sex::Female || sexes[b] == Sex::Female);
    };
    remove_edge_if(is_x, graph);

    auto vertex_range = boost::make_iterator_range(vertices(graph));
    for (auto v : vertex_range) {
        switch(sexes[v]) {
        case Sex::Female:
            clear_vertex(v, graph);
            ploidies[v] = 0;
            break;
        case Sex::Male:
            ploidies[v] = 1;
            break;
        default:
            if(out_degree(v, graph) != 0) {
                throw std::invalid_argument("Y-linked inheritance requires every individual to have a known sex.");
            }
        }
    }
}

void prune_xlinked(pedigree_graph::Graph &graph) {
    auto sexes  = get(boost::vertex_sex, graph);
    auto ploidies  = get(boost::vertex_ploidy, graph);

    auto is_y = [&](auto &&e) {
        if(!(get(boost::edge_type, graph, e) & mutk::detail::GERM_EDGE)) {
            return false;
        }
        auto a = source(e, graph);
        auto b = target(e, graph);
        return (sexes[a] == Sex::Male && sexes[b] == Sex::Male);
    };
    remove_edge_if(is_y, graph);

    auto vertex_range = boost::make_iterator_range(vertices(graph));
    for (auto v : vertex_range) {
        switch(sexes[v]) {
        case Sex::Female:
             break;
        case Sex::Male:
            ploidies[v] = 1;
            break;
        default:
            if(out_degree(v, graph) != 0) {
                throw std::invalid_argument("X-linked inheritance requires every individual to have a known sex.");
            }
        }
    }
}

void prune_wlinked(pedigree_graph::Graph &graph) {
    auto sexes  = get(boost::vertex_sex, graph);
    auto ploidies  = get(boost::vertex_ploidy, graph);

    auto is_z = [&](auto &&e) {
        if(!(get(boost::edge_type, graph, e) & mutk::detail::GERM_EDGE)) {
            return false;
        }
        auto a = source(e, graph);
        auto b = target(e, graph);
        return (sexes[a] == Sex::Male || sexes[b] == Sex::Male);
    };
    remove_edge_if(is_z, graph);

    auto vertex_range = boost::make_iterator_range(vertices(graph));
    for (auto v : vertex_range) {
        switch(sexes[v]) {
        case Sex::Male:
            clear_vertex(v, graph);
            ploidies[v] = 0;
            break;
        case Sex::Female:
            ploidies[v] = 1;
            break;
        default:
            if(out_degree(v, graph) != 0) {
                throw std::invalid_argument("W-linked inheritance requires every individual to have a known sex.");
            }
        }
    }
}

void prune_zlinked(pedigree_graph::Graph &graph) {
    auto sexes  = get(boost::vertex_sex, graph);
    auto ploidies  = get(boost::vertex_ploidy, graph);

    auto is_w = [&](auto &&e) {
        if(!(get(boost::edge_type, graph, e) & mutk::detail::GERM_EDGE)) {
            return false;
        }
        auto a = source(e, graph);
        auto b = target(e, graph);
        return (sexes[a] == Sex::Female && sexes[b] == Sex::Female);
    };
    remove_edge_if(is_w, graph);

    auto vertex_range = boost::make_iterator_range(vertices(graph));
    for (auto v : vertex_range) {
        switch(sexes[v]) {
        case Sex::Male:
             break;
        case Sex::Female:
            ploidies[v] = 1;
            break;
        default:
            if(out_degree(v, graph) != 0) {
                throw std::invalid_argument("Z-linked inheritance requires every individual to have a known sex.");
            }
        }
    }
}

void prune_maternal(pedigree_graph::Graph &graph) {
    auto sexes  = get(boost::vertex_sex, graph);
    auto ploidies  = get(boost::vertex_ploidy, graph);

    auto is_p = [&](auto &&e) {
        if(!(get(boost::edge_type, graph, e) & mutk::detail::GERM_EDGE)) {
            return false;
        }
        auto a = source(e, graph);
        return sexes[a] == Sex::Male;
    };
    remove_edge_if(is_p, graph);

    auto vertex_range = boost::make_iterator_range(vertices(graph));
    for (auto v : vertex_range) {
        ploidies[v] = 1;
    }
}

void prune_paternal(pedigree_graph::Graph &graph) {
    auto sexes  = get(boost::vertex_sex, graph);
    auto ploidies  = get(boost::vertex_ploidy, graph);

    auto is_m = [&](auto &&e) {
        if(!(get(boost::edge_type, graph, e) & mutk::detail::GERM_EDGE)) {
            return false;
        }
        auto a = source(e, graph);
        return sexes[a] == Sex::Male;
    };
    remove_edge_if(is_m, graph);

    auto vertex_range = boost::make_iterator_range(vertices(graph));
    for (auto v : vertex_range) {
        ploidies[v] = 1;
    }
}

void add_edges_to_pedigree_graph(const Pedigree &pedigree, pedigree_graph::Graph &graph) {
    auto has_tag = [](const Pedigree::Member &member, const std::string &tag) {
        for(auto &&a : member.tags) {
            if(boost::algorithm::iequals(a, tag)) {
                return true;
            }
        }
        return false;
    };

    // Add edges to graph
    auto vertex_range = boost::make_iterator_range(vertices(graph));
    for(auto v : vertex_range) {
        auto & member = pedigree.GetMember(v);
        if( has_tag(member, "founder") || ( !member.dad && !member.mom ) ) {
            continue;
        }
        auto ploidy = get(boost::vertex_ploidy, graph, v);

        if( ploidy == 0 /* clone */ ) {
            if(member.dad && member.mom) {
                throw std::invalid_argument("Unable to construct graph for pedigree; clone '"
                    + member.name + "' has two parents instead of one.");
            }
            
            // Determine which parent is the clone parent
            pedigree_graph::vertex_t parent;
            float len;
            if(member.dad) {
                parent = pedigree.LookupMemberPosition(*member.dad);
                len = member.dad_length.value_or(1.0);
            } else {
                parent = pedigree.LookupMemberPosition(*member.mom);
                len = member.mom_length.value_or(1.0);
            }
            if(parent >= pedigree.NumberOfMembers()) {
                throw std::invalid_argument("Unable to construct graph for pedigree; the clone parent of '" +
                    member.name + "' is unknown.");
            }
            // construct edge
            add_edge(parent, v, {len, mutk::detail::GERM_EDGE}, graph);

            // copy properties from parent to clone
            ploidy = get(boost::vertex_ploidy, graph, parent);
            put(boost::vertex_ploidy, graph, v, ploidy);
            auto sex = get(boost::vertex_sex, graph, parent);
            put(boost::vertex_sex, graph, v, sex);
        } else if(ploidy == 1 /* haploid/gamete */ ) {
            if(member.dad && member.mom) {
                throw std::invalid_argument("Unable to construct graph for pedigree; gamete '"
                    + member.name + "' has two parents instead of one.");
            }
            // Determine the parent
            pedigree_graph::vertex_t parent;
            float len;
            if(member.dad) {
                parent = pedigree.LookupMemberPosition(*member.dad);
                len = member.dad_length.value_or(1.0);
                if (get(boost::vertex_sex, graph, parent) == Sex::Female) {
                    throw std::invalid_argument("Unable to construct graph for pedigree; the father of '" +
                        member.name + "' is female.");
                }
            } else {
                parent = pedigree.LookupMemberPosition(*member.mom);
                len = member.mom_length.value_or(1.0);
                if (get(boost::vertex_sex, graph, parent) == Sex::Male) {
                    throw std::invalid_argument("Unable to construct graph for pedigree; the mother of '" +
                        member.name + "' is male.");
                }
            }
            if(parent >= pedigree.NumberOfMembers()) {
                throw std::invalid_argument("Unable to construct graph for pedigree; the parent of '" +
                    member.name + "' is unknown.");
            }
            // construct edge
            add_edge(parent, v, {len, mutk::detail::GERM_EDGE}, graph);
        } else {
            // Check for valid ids
            if(!member.dad) {
                throw std::invalid_argument("Unable to construct graph for pedigree; the father of '" +
                    member.name + "' is unspecified.");
            }
            if(!member.mom) {
                throw std::invalid_argument("Unable to construct graph for pedigree; the mother of '" +
                    member.name + "' is unspecified.");
            }
            pedigree_graph::vertex_t dad = pedigree.LookupMemberPosition(*member.dad);
            pedigree_graph::vertex_t mom = pedigree.LookupMemberPosition(*member.mom);
            float dad_len = member.dad_length.value_or(1.0);
            float mom_len = member.mom_length.value_or(1.0);

            if(dad >= pedigree.NumberOfMembers()) {
                throw std::invalid_argument("Unable to construct graph for pedigree; the father of '" +
                    member.name + "' is unknown.");
            }
            if(mom >= pedigree.NumberOfMembers()) {
                throw std::invalid_argument("Unable to construct graph for pedigree; the mother of '" +
                    member.name + "' is unknown.");
            }
            // Check for sanity
            if (get(boost::vertex_sex, graph, dad) == Sex::Female) {
                throw std::invalid_argument("Unable to construct graph for pedigree; the father of '" +
                    member.name + "' is female.");
            }
            if (get(boost::vertex_sex, graph, mom) == Sex::Male ) {
                throw std::invalid_argument("Unable to construct graph for pedigree; the mother of '" +
                    member.name + "' is male.");
            }
            add_edge(dad, v, {dad_len, mutk::detail::GERM_EDGE}, graph);
            add_edge(mom, v, {mom_len, mutk::detail::GERM_EDGE}, graph);
        }
    }
}

void construct_pedigree_graph(pedigree_graph::Graph &graph,
    const Pedigree &pedigree, const samples_t& known_samples,
    bool normalize_somatic_trees) {

    using namespace std;

    auto get_ploidy = [](const Pedigree::Member &member) -> int {
        using boost::algorithm::iequals;
        for(auto &&a : member.tags) {
            if(iequals(a, "haploid") || iequals(a, "gamete")
                || iequals(a, "p=1") || iequals(a, "ploidy=1")) {
                return 1;
            }
            if(iequals(a, "diploid") || iequals(a, "p=2") || iequals(a, "ploidy=2")) {
                return 2;
            }
        }
        // if a clone has no specified ploidy, set its ploidy to 0 and resolve it later
        for(auto &&a : member.tags) {
            if(iequals(a, "clone")) {
                return 0;
            }
        }
        // default to ploidy of 2
        return 2;
    };

    graph.clear();

    for(unsigned int j = 0; j < pedigree.NumberOfMembers(); ++j) {
        auto & member = pedigree.GetMember(j);
        auto v = add_vertex({member.name, {member.sex, {get_ploidy(member), VertexType::Germline}}}, graph);
        assert(v == j);
    }
    add_edges_to_pedigree_graph(pedigree, graph);

    // Add somatic branches and nodes
    for(int i=0;i<pedigree.NumberOfMembers();++i) {
        for(auto && sample : pedigree.GetMember(i).samples) {
            if(!mutk::detail::parse_newick(sample, graph, i,
                normalize_somatic_trees) ) {
                    throw std::invalid_argument("Unable to parse somatic data for individual '" +
                        pedigree.GetMember(i).name + "'.");
            }            
        }
    }

    // Identify and mark somatic nodes that connect to known samples
    auto vertex_range = boost::make_iterator_range(vertices(graph));
    auto types = get(boost::vertex_type, graph);
    auto labels = get(boost::vertex_label, graph);
    std::unordered_set<std::string> known(known_samples.begin(), known_samples.end());
    for(auto v : vertex_range) {
        if(types[v] != VertexType::Somatic) {
            continue;
        }
        if(known.find(labels[v]) != known.end()) {
            types[v] = VertexType::Sample;
        }
    }
}

void update_edge_lengths(pedigree_graph::Graph &graph,
        double mu_germ, double mu_soma) {
    auto edge_types = get(boost::edge_type, graph);
    auto lengths = get(boost::edge_length, graph);

    auto range = boost::make_iterator_range(edges(graph));
    for (auto && e : range) {
        if(edge_types[e] & mutk::detail::GERM_EDGE) {
            lengths[e] *= mu_germ;
        } else {
            lengths[e] *= mu_soma;
        }
    }
}

void simplify(pedigree_graph::Graph &graph) {
    using boost::make_iterator_range;

    auto edge_types = get(boost::edge_type, graph);
    auto lengths = get(boost::edge_length, graph);
    auto types = get(boost::vertex_type, graph);
    auto ploidies = get(boost::vertex_ploidy, graph);

    // topologically sort members in reverse order
    std::vector<pedigree_graph::vertex_t> rev_topo_order;
    topological_sort(graph, std::back_inserter(rev_topo_order));
    auto topo_order = boost::adaptors::reverse(rev_topo_order);

    // Clear all leaf vertexes that are not samples, starting from the tips
    for(auto v : rev_topo_order) {
        if(out_degree(v, graph) == 0 && types[v] != VertexType::Sample) {
            clear_vertex(v, graph);
        }
    }
    // unlink founder nodes that are summed out anyways
    for(auto && v : topo_order) {
        if(types[v] != VertexType::Germline) {
            continue;
        }
        auto edge_range = boost::make_iterator_range(in_edges(v,graph));
        if(edge_range.empty()) {
            continue;
        }
        bool all = boost::algorithm::all_of(edge_range, [&](auto &&e){
            auto p = source(e, graph);
            return (degree(p, graph) == 1);
        });
        if (all) {
            clear_in_edges(v, graph);
        }
    }

    // bypass nodes that have one out_edge and their descendants
    for(auto && v : topo_order) {
        if(in_degree(v,graph) == 0 || out_degree(v,graph) != 1) {
            continue;
        }
        auto in_edge_range = boost::make_iterator_range(in_edges(v,graph));
        auto out_edge = *out_edges(v,graph).first;
        auto child = target(out_edge, graph);
        if(in_degree(child,graph)+in_degree(v,graph)-1 > 2) {
            continue;
        }
        if(ploidies[child] != ploidies[v]) {
            continue;
        }
        auto len = lengths[out_edge];
        EdgeType otype = edge_types[out_edge];
        for(auto &&e : in_edge_range) {
            EdgeType etype = static_cast<EdgeType>( otype | edge_types[e]);
            auto grand = source(e, graph);
            add_edge(grand, child, {len+lengths[e], etype}, graph);
        }
        clear_vertex(v, graph);
    }
}

pedigree_graph::Graph finalize(const pedigree_graph::Graph &input) {
    // Construct new graph
    pedigree_graph::Graph output;
    auto types = get(boost::vertex_type, input);
    // topologically sort members
    std::vector<pedigree_graph::vertex_t> topo_order, vertex_order;
    topological_sort(input, std::back_inserter(topo_order));
    // Founders
    std::copy_if(topo_order.rbegin(), topo_order.rend(), std::back_inserter(vertex_order),
        [&](auto v) {
            return in_degree(v, input) == 0 &&
                out_degree(v, input) > 0 &&
                types[v] == VertexType::Germline;
        });
    // Germline
    std::copy_if(topo_order.rbegin(), topo_order.rend(), std::back_inserter(vertex_order),
        [&](auto v) {
            return in_degree(v, input) > 0 && types[v] == VertexType::Germline;
        });
    // Somatic
    std::copy_if(topo_order.rbegin(), topo_order.rend(),
        std::back_inserter(vertex_order), [&](auto v) {
            return degree(v, input) > 0 && types[v] == VertexType::Somatic;
        });
    // Samples
    std::copy_if(topo_order.rbegin(), topo_order.rend(),
        std::back_inserter(vertex_order), [&](auto v) {
            return degree(v, input) > 0 && types[v] == VertexType::Sample;
        });

    std::vector<pedigree_graph::vertex_t> map_in_to_out(num_vertices(input),-1);

    auto local_vertex_map = get(boost::vertex_all, input);
    auto graph_vertex_map = get(boost::vertex_all, output);
    for(auto &&v : vertex_order) {
        auto w = add_vertex(output);
        map_in_to_out[v] = w;
        put(graph_vertex_map, w, get(local_vertex_map, v));
    }

    auto input_edge_map = get(boost::edge_all, input);
    auto output_edge_map = get(boost::edge_all, output); 
    auto edge_range = boost::make_iterator_range(edges(input));
    for(auto &&e : edge_range) {
        auto new_src = map_in_to_out[source(e, input)];
        auto new_tgt = map_in_to_out[target(e, input)];
        assert(new_src != -1 && new_tgt != -1);
        auto [f,b] = add_edge(new_src, new_tgt, output);
        put(output_edge_map, f, get(input_edge_map, e));
    }

    auto olabels = get(boost::vertex_label, output);
    auto otypes = get(boost::vertex_type, output);
    for(auto &&v : boost::make_iterator_range(vertices(output))) {
        if(otypes[v] == VertexType::Germline) {
            olabels[v] += "/z";
        } else if(otypes[v] == VertexType::Somatic) {
            olabels[v] += "/t";
        }
    }

    return output;
}

std::vector<mutk::potential_t>
create_potentials(const pedigree_graph::Graph &graph) {
    using PotentialType = mutk::PotentialType;

    auto ploidies = get(boost::vertex_ploidy, graph);
    auto lengths = get(boost::edge_length, graph);

    std::vector<mutk::potential_t> ret;
    
    auto vertex_range = boost::make_iterator_range(vertices(graph));   
    // add samples
    for(auto v : vertex_range) {
        // sanity checks
        assert(ploidies[v] == 1 || ploidies[v] == 2);
        assert(in_degree(v,graph) <= 2);
        if(out_degree(v,graph) == 0) {
            auto pot_type = (ploidies[v] == 1) ? PotentialType::LikelihoodHaploid
                                               : PotentialType::LikelihoodDiploid;
            ret.emplace_back(pot_type, IndyId(v));
        }
    }
    // add founders
    for(auto v : vertex_range) {
        if(in_degree(v,graph) == 0) {
            auto pot_type = (ploidies[v] == 1) ? PotentialType::FounderHaploid
                                               : PotentialType::FounderDiploid;
            ret.emplace_back(pot_type, IndyId(v));
            // founders can't also be samples, as it ruins assumptions.
            assert(out_degree(v,graph) > 0);
        }
    }
    // add transitions
    for(auto v : vertex_range) {
        if(in_degree(v,graph) == 1) {
            auto [ei, ee] = in_edges(v,graph);
            auto  par1 = source(*ei,graph);
            float len1 = lengths[*ei];
            // identify the type of transition
            PotentialType pot_type;
            if(ploidies[par1] == 2) {
                // diploid -> diploid (clone) or
                // diploid -> haploid (gamete)
                pot_type = (ploidies[v] == 1) ? PotentialType::GameteDiploid
                                              : PotentialType::CloneDiploid;
            } else {
                // haploid -> haploid (clone)
                assert(ploidies[v] == 1);
                pot_type = PotentialType::CloneHaploid;
            }
            ret.emplace_back(pot_type, IndyId(v), IndyId(par1), len1);
        } else if(in_degree(v,graph) == 2) {
            auto [ei, ee] = in_edges(v,graph);
            auto  par1 = source(*ei,graph);
            float len1 = lengths[*ei];
            ++ei;
            auto  par2 = source(*ei,graph);
            float len2 = lengths[*ei];
            // identify the type of transition
            PotentialType pot_type;
            if(par1 == par2) {
                // selfing
                pot_type = (ploidies[v] == 1) ? PotentialType::ChildSelfingHaploid
                                              : PotentialType::ChildSelfingDiploid;
            } else if(ploidies[par1] == 2) {
                pot_type = (ploidies[par2] == 1) ? PotentialType::ChildDiploidHaploid
                                                 : PotentialType::ChildDiploidDiploid;
            } else {
                pot_type = (ploidies[par2] == 1) ? PotentialType::ChildHaploidHaploid
                                                 : PotentialType::ChildHaploidDiploid;                
            }
            ret.emplace_back(pot_type, IndyId(v), IndyId(par1), len1, IndyId(par2), len2);
        }
    }
    return ret;
}

using record_t = std::pair<int,pedigree_graph::vertex_t>;

struct record_cmp_t {
    bool operator()(const record_t &a, const record_t &b) const {
        if(a.first == b.first) {
            return a.second < b.second;
        } else {
            return a.first > b.first;
        }
    }
};

using heap_t = boost::heap::d_ary_heap<record_t,
    boost::heap::arity<2>, boost::heap::mutable_<true>,
    boost::heap::compare<record_cmp_t>>;

// Almond and Kong (1991) Optimality Issues in Constructing a Markov Tree from Graphical Models.
//     Research Report 329. University of Chicago, Dept. of Statistics
std::pair<std::vector<pedigree_graph::vertex_t>,
std::vector<neighbors_t>>
triangulate_graph(const pedigree_graph::Graph &graph) {
    using vertex_t = pedigree_graph::vertex_t;

    auto vertex_range = boost::make_iterator_range(vertices(graph));

    // Identify Closures
    std::vector<neighbors_t> neighbors(num_vertices(graph));
    for(auto v: vertex_range) {
        auto [it1, itend] = in_edges(v, graph);
        for(; it1 != itend; ++it1) {
            auto it2 = it1;
            // connect v to each of its parents
            auto p1 = source(*it1, graph);
            neighbors[v].insert(p1);
            neighbors[p1].insert(v);
            for(++it2; it2 != itend; ++it2) {
                // connect each of the parents to one another
                auto p2 = source(*it2, graph);
                neighbors[p1].insert(p2);
                neighbors[p2].insert(p1);
            }
        }
    }

    auto fill_in_count = [&](auto v) {
        int fill = 0;
        const auto& k = neighbors[v];
        for(auto it1 = k.begin(); it1 != k.end(); ++it1) {
            auto it2 = it1;
            for(++it2; it2 != k.end(); ++it2) {
                if(!neighbors[*it1].contains(*it2)) {
                    fill += 1;
                }
                //sanity check
                assert(neighbors[*it1].contains(*it2) == neighbors[*it2].contains(*it1));
            }
        }
        return fill;
    };

    heap_t priority_queue;

    // create priority queue based on fill-in
    std::vector<heap_t::handle_type> handles(num_vertices(graph));
    for(auto v : boost::make_iterator_range(vertices(graph))) {
        int f = fill_in_count(v);
        handles[v] = priority_queue.push({f, v});
    }

    std::vector<vertex_t> elim_order;

    while(!priority_queue.empty()) {
        // chose next vertex to eliminate and remove it from the queue
        auto [score, v] = priority_queue.top();
        priority_queue.pop();
        handles[v] = heap_t::handle_type{};

        // record the vertex
        elim_order.push_back(v);

        // Update cliques
        auto &k = neighbors[v];
        for(auto &a : k) {
            if(score > 0) {
                for(auto &b : k) {
                    neighbors[a].insert(b);
                }                
            }
            neighbors[a].erase(v);
        }
        // Update the priority queue
        for(auto &a : k) {
            (*handles[a]).first = fill_in_count(a);
            priority_queue.update(handles[a]);
        }
    }
    // DEBUG: print elimination information
    // for(auto &&v : elim_order) {
    //     std::cerr << "Eliminate vertex " << get(boost::vertex_label, graph, v)
    //         << " clique is { ";
    //     std::cerr << get(boost::vertex_label, graph, v);
    //     for(auto &&a : neighbors[v]) {
    //          std::cerr << ", " << get(boost::vertex_label, graph, a);
    //     }
    //     std::cerr << " }\n";
    // }

    return {elim_order,neighbors};
}

junction_tree::Graph create_junction_tree(
    const std::vector<pedigree_graph::vertex_t> &elim_order,
    const std::vector<neighbors_t> &neighbors) {

    junction_tree::Graph ret;
    auto node_labels = get(boost::vertex_label, ret);

    for(auto v : boost::adaptors::reverse(elim_order)) {
        // add vertex to junction_tree
        junction_tree::vertex_t j;

        auto node_range = boost::make_iterator_range(vertices(ret));
        
        if(neighbors[v].empty()) {
            j = add_vertex({{v}}, ret);
            continue;
        }
        auto nit = boost::range::find_if(node_range, [&](const auto &n) {
                return (neighbors[v] == node_labels[n]);
        });
        if(nit != boost::end(node_range)) {
            // neighbors of v already exists as a vertex
            // mark as intersection node and create edge
            j = add_vertex(ret);
            add_edge(*nit, j, ret);
        } else {
            using boost::adaptors::filtered;
            // find the smallest vertex that contains neighbors as a subset.
            auto is_subset = [&] (const auto &x) {
                return boost::range::includes(node_labels[x], neighbors[v]);
            };
            auto cmp_size = [&] (const auto &a, const auto &b) {
                return node_labels[a].size() < node_labels[b].size();
            };
            nit = boost::range::min_element(node_range | filtered(is_subset), cmp_size).base();
            assert(nit != boost::end(node_range));
            // create an intersection node for the neighbors
            auto p = add_vertex(ret);
            node_labels[p] = neighbors[v];
            add_edge(*nit, p, ret);
            // create a new node
            j = add_vertex(ret);
            add_edge(p, j, ret);
        }
        // update vertex properties
        node_labels[j] = neighbors[v];
        node_labels[j].insert(v);

        // attach singleton vertex
        auto k = add_vertex({{v}}, ret);
        add_edge(j, k, ret);
    }

    // Topologically Sort Junction Tree
    std::vector<junction_tree::vertex_t> ordered_junction_vertices;
    boost::topo_sort_visitor sort_visitor(std::back_inserter(ordered_junction_vertices));
    boost::undirected_dfs(ret, boost::root_vertex(junction_tree::vertex_t(0))
        .visitor(sort_visitor)
        .edge_color_map(get(boost::edge_color, ret))
    );
    return ret;
}
} // anon namespace 


mutk::GeneralPeelingVertex::data_t
mutk::GeneralPeelingVertex::MakeMetadata(const std::vector<PeelingVertex::label_t> &labels,
    std::size_t buffer_index) const {
    data_t ret;
    ret.index = buffer_index;
    // fill with zeros.
    ret.shape.assign(labels_.size(), 0);

    // Go through the labels for this vertex and see if they are in the Metadata labels
    // ret.shape[i] is 1 if the label of axis i is in the labels of this vertex
    for(size_t i = 0; i < labels_.size(); ++i) {
        auto it = boost::range::find(labels, labels_[i]);
        if(it != labels.end()) {
            // our metadata contains the label
            ret.shape[i] = 1;
        }
    }
    return ret;
}

// LCOV_EXCL_START
namespace {
class Test_MakeMetadata : public mutk::GeneralPeelingVertex {
    using mutk::GeneralPeelingVertex::GeneralPeelingVertex;
TEST_CASE_CLASS("GeneralPeelingVertex-MakeMetadata") {
    Test_MakeMetadata v({1,0});

    using shape_t = mutk::shape_t;

    auto d = v.MakeMetadata({0}, 0);
    CHECK(d.index == 0);
    CHECK(d.shape == shape_t{0,1});

    v.AddLocal({1,0}, 0);
    v.AddInput({0},  1);
    v.AddOutput({1}, 2);

    CHECK(v.local_data_.index == 0);
    CHECK(v.local_data_.shape == shape_t{1,1});

    CHECK(v.input_data_.size() == 1);
    CHECK(v.input_data_[0].index == 1);
    CHECK(v.input_data_[0].shape == shape_t{0,1});

    CHECK(v.output_data_.index == 2);
    CHECK(v.output_data_.shape == shape_t{1});
}};}
// LCOV_EXCL_STOP

void mutk::GeneralPeelingVertex::AddLocal(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) {
    local_data_ = MakeMetadata(labels, buffer_index);
}

void mutk::GeneralPeelingVertex::AddOutput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) {
    output_data_ = MakeMetadata(labels, buffer_index);
    // for output_data, convert shape to a vector of the axes that need to be 
    // summed out when creating the output
    shape_t shape;
    for(size_t i=0; i < output_data_.shape.size(); ++i) {
        // if axis is missing in the output, add it to the list that needs to be removed.
        if(output_data_.shape[i] == 0) {
            shape.push_back(i);
        }
    }
    // update shape
    output_data_.shape = shape;
}

void mutk::GeneralPeelingVertex::AddInput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) {
    input_data_.push_back(MakeMetadata(labels, buffer_index));
}

namespace {
    mutk::shape_t calc_reshape(mutk::shape_t a, const mutk::shape_t& b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), [&](auto x, auto y) {
        return (y == 0) ? 1 : x;
    });
    return a;
}
} // anon namespace

// Forward Algorithm:
void mutk::GeneralPeelingVertex::Forward(PeelingVertex::workspace_t &work) const {
    // copy local data to local buffer
    temporary_ = work.stack[local_data_.index];
    for(int i=0;i<input_data_.size(); ++i) {
        auto dims = calc_reshape(temporary_.shape(), input_data_[i].shape);
        temporary_ *= xt::reshape_view(work.stack[input_data_[i].index], dims);
    }
    work.stack[output_data_.index] = xt::sum(temporary_, output_data_.shape);
}

// LCOV_EXCL_START
TEST_CASE("GeneralPeelingVertex-Forward") {
    using tensor_t = mutk::tensor_t;
    mutk::RelationshipGraph::workspace_t work;
    work.scale = 0.0;
    work.stack.resize(10);

    SUBCASE("1 Dimensional Vertex") {
        // Local
        work.stack[0] = tensor_t({3}, 1.0);
        // Input
        work.stack[1] = tensor_t({3}, 0.5);
        // Output
        work.stack[2].resize({});

        tensor_t expected = work.stack[0] * work.stack[1];

        mutk::GeneralPeelingVertex v({0});
        v.AddLocal({0}, 0);
        v.AddInput({0}, 1);
        v.AddOutput({0},2);
        v.Forward(work);

        REQUIRE(work.stack[2].size() == 3);
        for(int k=0; k < work.stack[2].size(); ++k) {
            CAPTURE(k);
            CHECK(work.stack[2][k] == expected[k]);            
        }
    }
    SUBCASE("2 Dimensional Vertex") {
        // Local
        work.stack[0] = {{0.8,0.2,0.3},
                         {0.1,0.6,0.3},
                         {0.1,0.2,0.4}};
        // Input
        work.stack[1] = {0.0,0.1,0.9};
        // Output
        work.stack[2].resize({});

        // Expected output
        tensor_t expected = xt::linalg::dot(work.stack[0], work.stack[1]);
        // expected = 0.29,0.33,0.38

        mutk::GeneralPeelingVertex v({0,1});
        v.AddLocal({0,1}, 0);
        v.AddInput({1},  1);
        v.AddOutput({0}, 2);
        v.Forward(work);

        REQUIRE(work.stack[2].size() == 3);
        for(int k=0; k < work.stack[2].size(); ++k) {
            CAPTURE(k);
            CHECK(work.stack[2](k) == expected(k));            
        }
    }
}
// LCOV_EXCL_STOP