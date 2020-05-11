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

#include <mutk/relationship_graph.hpp>
#include <mutk/detail/graph.hpp>

#include <queue>
#include <algorithm>
#include <unordered_set>

#include <boost/range/algorithm/find.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/container/flat_set.hpp>

#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/heap/d_ary_heap.hpp>

namespace {

using Sex = mutk::Pedigree::Sex;

using samples_t = mutk::RelationshipGraph::samples_t;

using InheritanceModel = mutk::InheritanceModel;
using Pedigree = mutk::Pedigree;

namespace pedigree_graph = mutk::detail::pedigree_graph;
using VertexType = mutk::detail::VertexType;
using EdgeType = mutk::detail::EdgeType;

void construct_pedigree_graph(pedigree_graph::Graph &graph,
    const Pedigree &pedigree, const samples_t& known_samples,
    bool normalize_somatic_trees);

void update_edge_lengths(pedigree_graph::Graph &graph, double mu_meiotic, double mu_somatic);
void simplify(pedigree_graph::Graph &graph);
void prune(pedigree_graph::Graph &graph, InheritanceModel model);

pedigree_graph::Graph finalize(const pedigree_graph::Graph &input);

void peeling_order(const pedigree_graph::Graph &graph);

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

bool mutk::RelationshipGraph::Construct(const Pedigree& pedigree,
        const samples_t& known_samples, InheritanceModel model,
        double mu, double mu_somatic,
        bool normalize_somatic_trees) {
    using namespace std;

    inheritance_model_ = model;

    // Construct a boost::graph of the pedigree and somatic information
    pedigree_graph::Graph graph;

    for(auto v : known_samples) {
        std::cout << v << "\n";
    }

    construct_pedigree_graph(graph, pedigree, known_samples, normalize_somatic_trees);
 
    // Multiply edge lengths by mutation rates
    update_edge_lengths(graph, mu, mu_somatic);

    pedigree_graph::print_graph(graph);

    // Remove edges that are non-informative
    simplify(graph);

    // Prune pedigree
    prune(graph, inheritance_model_);

    pedigree_graph::print_graph(graph);

    // Sort and eliminate cleared vertices
    graph_ = finalize(graph);
    
    return true;
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

    auto print_vertex = [&](auto &&v, const char *str) {
        os << labels[v] << "\t"
           << str << "\t"
           << sex(sexes[v]) << "\t"
           << ploidies[v] << "\t";
        auto [ei,ee] = in_edges(v, graph_);
        for(auto it = ei; it != ee; ++it) {
            if(it != ei) {
                os << ";";
            }
            os << labels[source(*it,graph_)] << ":"
                << get(boost::edge_length, graph_, *it);
        }
        os << "\n";
    };

    os << "Name\tType\tSex\tPloidy\tSource\n";

    for(auto v : vertex_range) {
        if(in_degree(v,graph_) == 0) {
            // add samples
            print_vertex(v, "FOUNDER");
        }
    }
    for(auto v : vertex_range) {
        if(in_degree(v,graph_) > 0 && types[v] == VertexType::Germline) {
            // add samples
            print_vertex(v, "GERMLINE");
        }
    }
    for(auto v : vertex_range) {
        if(in_degree(v,graph_) > 0 && types[v] == VertexType::Somatic) {
            // add samples
            print_vertex(v, "SOMA");
        }
    }
    for(auto v : vertex_range) {
        if(in_degree(v,graph_) > 0 && types[v] == VertexType::Sample) {
            // add samples
            print_vertex(v, "SAMPLE");
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
            // Selfing is not supported
            if (dad == mom) {
                throw std::invalid_argument("Unable to construct graph for pedigree; selfing is not supported; "
                    "father and mother of '" + member.name + "' are the same.");
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
    using Sex = mutk::Pedigree::Sex;

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

    // bypass nodes that have one out_edge and their descendents
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
        std::cout << "vertex " << v << " maps to vertex " << w << "\n";
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

    return output;
}

// Almond and Kong (1991) Optimality Issues in Constructing a Markov Tree from Graphical Models.
//     Research Report 329. University of Chicago, Dept. of Statistics

template<int N_>
//using small_vector_t = boost::container::small_vector<pedigree_graph::vertex_t, N_>;
using small_vector_t = std::vector<pedigree_graph::vertex_t>;

template<int N_>
using flat_set_t = boost::container::flat_set<pedigree_graph::vertex_t,
        std::less<pedigree_graph::vertex_t>,
        small_vector_t<N_>>;

void peeling_order(const pedigree_graph::Graph &graph) {
    using vertex_t = pedigree_graph::vertex_t;

    auto types = get(boost::vertex_type, graph);

    auto vertex_range = boost::make_iterator_range(vertices(graph));

    // Identify vertex dependencies
    //   - dependencies[i] contains the in-vertexes of i.
    using depends_t = small_vector_t<2>;
    std::vector<depends_t> depends(num_vertices(graph));
    for(auto v : vertex_range) {
        auto in_edge_range = boost::make_iterator_range(in_edges(v,graph));
        for(auto e : in_edge_range) {
            depends[v].push_back(source(e,graph));
        }
        boost::sort(depends[v]);
    }
    // Factorize the probability distribution into potentials
    using potential_t = small_vector_t<4>;
    std::vector<potential_t> potentials(num_vertices(graph));
    for(auto v : vertex_range) {
        if(out_degree(v,graph) == 0) {
            // add samples
            potentials.push_back({v});
        }
        if(in_degree(v,graph) == 0) {
            // Founders
            potentials.push_back({v});
        } else {
            // Non-Founders
            potentials.push_back({v});
            potentials.back().insert(potentials.back().end(),
                depends[v].begin(), depends[v].end());

        }
    }

    // Identify Cliques
    using clique_t = flat_set_t<8>;
    std::vector<clique_t> cliques(num_vertices(graph));
    for(auto && p : potentials) {
        for(auto && v : p) {
            for(auto && w : p) {
                cliques[v].insert(w);
            }
        }
    }
    std::cout << "=== CLIQUES ===\n";
    for(auto &k : cliques) {
        for(auto &v : k) {
            std::cout << v << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";

    auto fill_in_count = [&](auto v) {
        int fill = 0;
        const auto& k = cliques[v];
        for(auto it1 = k.begin(); it1 != k.end(); ++it1) {
            auto it2 = it1;
            for(++it2; it2 != k.end(); ++it2) {
                if(!cliques[*it1].contains(*it2)) {
                    fill += 1;
                }
                //sanity check
                assert(cliques[*it1].contains(*it2) == cliques[*it2].contains(*it1));
            }
        }
        return fill;
    };

    using record_t = std::pair<int,vertex_t>;
    using heap_t = boost::heap::d_ary_heap<record_t,
        boost::heap::arity<2>, boost::heap::mutable_<true>,
        boost::heap::compare<std::greater<record_t>>>;

    heap_t priority_queue;

    // create priority queue based on fill-in
    std::vector<heap_t::handle_type> handles(num_vertices(graph));
    for(auto v : boost::make_iterator_range(vertices(graph))) {
        int f = fill_in_count(v);
        handles[v] = priority_queue.push({f, v});
    }

    while(!priority_queue.empty()) {
        for(auto it = priority_queue.ordered_begin();
            it != priority_queue.ordered_end(); ++it) {
            std::cout << it->second << " " << it->first << "\n";
        }

        auto [score, v] = priority_queue.top();
        priority_queue.pop();
        handles[v] = heap_t::handle_type{};
        std::cout << "Eliminating Vertex " << v << "\n";

        auto &k = cliques[v];

        std::cout << "    Clique: ";
        for(auto &a : k) {
            std::cout << " " << a;
        }
        std::cout << "\n\n";

        // fill in cliques
        for(auto &a : k) {
            if(a == v) {
                continue;
            }
            if(score > 0) {
                for(auto &b : k) {
                    cliques[a].insert(b);
                }                
            }
            cliques[a].erase(v);
        }
        for(auto &a : k) {
            if(a == v) {
                continue;
            }
            (*handles[a]).first = fill_in_count(a);
            priority_queue.update(handles[a]);
        }
    }
}

} // namespace 

