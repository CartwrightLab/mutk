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
#include <boost/algorithm/cxx11/all_of.hpp>

#include <boost/graph/topological_sort.hpp>

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

    construct_pedigree_graph(graph, pedigree, known_samples, normalize_somatic_trees);

    std::cout << "====Pedigree Graph====\n";
    pedigree_graph::print_graph(graph);
 
    // Multiply edge lengths by mutation rates
    update_edge_lengths(graph, mu, mu_somatic);

    // Remove edges that are non-informative
    simplify(graph);

    std::cout << "====Pedigree Graph====\n";
    pedigree_graph::print_graph(graph);

    // Prune pedigree
    prune(graph, inheritance_model_);

    std::cout << "====Pedigree Graph====\n";
    pedigree_graph::print_graph(graph);
    
    return true;
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

    // Construct temporary graph
    pedigree_graph::Graph local_graph;

    for(unsigned int j = 0; j < pedigree.NumberOfMembers(); ++j) {
        auto & member = pedigree.GetMember(j);
        auto v = add_vertex({member.name, {member.sex, {get_ploidy(member), VertexType::Germline}}}, local_graph);
        assert(v == j);
    }
    add_edges_to_pedigree_graph(pedigree, local_graph);

    // topologically sort members
    std::vector<pedigree_graph::vertex_t> topo_order, vertex_order;

    topological_sort(local_graph, std::back_inserter(topo_order));

    // sort founders before everyone else
    std::copy_if(topo_order.rbegin(), topo_order.rend(), std::back_inserter(vertex_order),
        [&](auto v) {
            return in_degree(v, local_graph) == 0;
        });
    std::copy_if(topo_order.rbegin(), topo_order.rend(), std::back_inserter(vertex_order),
        [&](auto v) {
            return in_degree(v, local_graph) > 0;
        });

    // Copy local_graph to graph using vertex_order
    graph.clear();

    std::vector<pedigree_graph::vertex_t> map_local_to_out(num_vertices(local_graph));

    auto local_vertex_map = get(boost::vertex_all, local_graph);
    auto graph_vertex_map = get(boost::vertex_all, graph);
    for(auto &&v : vertex_order) {
        auto w = add_vertex(graph);
        map_local_to_out[v] = w;
        put(graph_vertex_map, w, get(local_vertex_map, v));
    }

    auto local_edge_map = get(boost::edge_all, local_graph);
    auto graph_edge_map = get(boost::edge_all, graph); 
    auto edge_range = boost::make_iterator_range(edges(local_graph));
    for(auto &&e : edge_range) {
        auto new_src = map_local_to_out[source(e, local_graph)];
        auto new_tgt = map_local_to_out[target(e, local_graph)];
        auto [f,b] = add_edge(new_src, new_tgt, graph);
        put(graph_edge_map, f, get(local_edge_map, e));
    }

    // Add somatic branches and nodes
    for(int i=0;i<pedigree.NumberOfMembers();++i) {
        for(auto && sample : pedigree.GetMember(i).samples) {
            if(!mutk::detail::parse_newick(sample, graph, map_local_to_out[i],
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

    // topologically sort members in reverse order
    std::vector<pedigree_graph::vertex_t> rev_topo_order;
    topological_sort(graph, std::back_inserter(rev_topo_order));
    auto topo_order = boost::adaptors::reverse(rev_topo_order);

    // Clear all leaf vertexes that are not samples, starting from the tips
    for(auto v : rev_topo_order) {
        if(out_degree(v, graph) == 0 && types[v] != VertexType::Sample) {
            std::cout << "clearing edges " << v << "\n";
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
            std::cout << "clearing in edges " << v << "\n";
            clear_in_edges(v, graph);
        }
    }

    // bypass nodes that have one out_edge
    for(auto && v : topo_order) {
        if(in_degree(v,graph) == 0 || out_degree(v,graph) != 1) {
            continue;
        }
        auto in_edge_range = boost::make_iterator_range(in_edges(v,graph));
        auto out_edge = *out_edges(v,graph).first;
        auto child = target(out_edge, graph);
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

} // namespace 

