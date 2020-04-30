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

#include <queue>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/adaptor/reversed.hpp>


// Install boost graph properties
namespace boost {
enum edge_length_t { edge_length };
enum edge_type_t { edge_type };
enum vertex_sex_t { vertex_sex };
enum vertex_ploidy_t {vertex_ploidy };
enum vertex_label_t { vertex_label };
enum vertex_library_label_t { vertex_library_label };
enum vertex_type_t { vertex_type };

enum edge_family_t { edge_family };
enum vertex_group_t { vertex_group };

BOOST_INSTALL_PROPERTY(edge, length);
BOOST_INSTALL_PROPERTY(edge, type);
BOOST_INSTALL_PROPERTY(vertex, sex);
BOOST_INSTALL_PROPERTY(vertex, ploidy);
BOOST_INSTALL_PROPERTY(vertex, label);
BOOST_INSTALL_PROPERTY(vertex, library_label);
BOOST_INSTALL_PROPERTY(vertex, type);

BOOST_INSTALL_PROPERTY(edge, family);
BOOST_INSTALL_PROPERTY(vertex, group);
}

namespace {

enum struct EdgeType : std::size_t {
    Spousal, Maternal, Paternal, Mitotic
};
enum struct VertexType : std::size_t {
    Germline, Somatic
};
using Sex = mutk::Pedigree::Sex;

using VertexGroupProp = boost::property<boost::vertex_group_t, std::size_t>;
using VertexLibraryLabelProp =  boost::property<boost::vertex_library_label_t, std::string, VertexGroupProp>;
using VertexPloidyProp = boost::property<boost::vertex_ploidy_t, int, VertexLibraryLabelProp>;
using VertexSexProp = boost::property<boost::vertex_sex_t, Sex, VertexPloidyProp>;
using VertexTypeProp = boost::property<boost::vertex_type_t, VertexType,VertexSexProp>;
using VertexLabelProp = boost::property<boost::vertex_label_t, std::string, VertexTypeProp>;

using EdgeFamilyProp = boost::property<boost::edge_family_t, std::size_t>;
using EdgeLengthProp = boost::property<boost::edge_length_t, float, EdgeFamilyProp>;
using EdgeTypeProp = boost::property<boost::edge_type_t, EdgeType, EdgeLengthProp>;

using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
        VertexLabelProp, EdgeTypeProp>;

using vertex_t = boost::graph_traits<Graph>::vertex_descriptor vertex_t;
using edge_t = boost::graph_traits<Graph>::edge_descriptor;

static_assert(std::is_integral<vertex_t>::value,
    "vertex_t is not an integral type, this violates many assumptions that have been made.");

using samples_t = mutk::RelationshipGraph::samples_t;

using Graph = mutk::detail::graph::Graph;
using vertex_t = mutk::detail::graph::vertex_t;
using InheritanceModel = mutk::InheritanceModel;
using Pedigree = mutk::Pedigree;

std::pair<vertex_t,vertex_t> parse_pedigree_table(Graph &pedigree_graph, const Pedigree &pedigree,
        bool normalize_somatic_trees);

bool parse_newick(const std::string &text, int root, Graph *graph, bool normalize);

void add_samples_to_graph(Graph &pedigree_graph, const libraries_t &libs);
void update_edge_lengths(Graph &pedigree_graph,
        double mu_meiotic, double mu_somatic, double mu_library);
void simplify_pedigree(Graph &pedigree_graph);

void prune_pedigree(Graph &pedigree_graph, InheritanceModel model);

void prefix_vertex_labels(Graph &pedigree_graph);

const std::pair<std::string, InheritanceModel> inheritance_keys[] = {
    {"", InheritanceModel::Invalid},
    {"AUTOSOMAL", InheritanceModel::Autosomal},
    {"MATERNAL", InheritanceModel::Maternal},
    {"PATERNAL", InheritanceModel::Paternal},
    {"X-LINKED", InheritanceModel::XLinked},
    {"Y-LINKED", InheritanceModel::YLinked},
    {"W-LINKED", InheritanceModel::WLinked},
    {"Z-LINKED", InheritanceModel::ZLinked},
    {"MITOCHONDRIAL", InheritanceModel::Maternal},
    {"XLINKED", InheritanceModel::XLinked},
    {"YLINKED", InheritanceModel::YLinked},
    {"WLINKED", InheritanceModel::WLinked},
    {"ZLINKED", InheritanceModel::ZLinked}
};

} // namespace

InheritanceModel mutk::inheritance_model(const std::string &pattern) {
    InheritanceModel model = utility::key_switch_tuple(pattern, inheritance_keys,
                                                inheritance_keys[0]).second;
    if (model == InheritanceModel::Invalid){
        throw std::invalid_argument("Inheritance model '" + pattern
            + "' is not supported. Supported values are: "
            "[autosomal, mitochondrial, paternal, x-linked, y-linked, w-linked, z-linked]");
    }
    return model;
}

std::string mutk::to_string(InheritanceModel model) {
    if(model != InheritanceModel::Invalid) {
        for(auto &&a : inheritance_keys) {
            if(std::get<1>(a) != model) {
                continue;
            }
            return std::get<0>(a);
        }
    }
    throw std::invalid_argument("Unable to convert model '" + std::to_string((int)model)
            + "' to string.");
    return {};   
}

bool mutk::RelationshipGraph::Construct(const Pedigree& pedigree,
        const samples_t& known_samples, InheritanceModel model,
        double mu, double mu_somatic, double mu_library,
        bool normalize_somatic_trees) {
    using namespace std;

    inheritance_model_ = model;

    // Construct a boost::graph of the pedigree and somatic information
    Graph pedigree_graph;

    first_founder_ = 0;
    std::tie(first_nonfounder_,first_somatic_) = parse_pedigree_table(pedigree_graph, pedigree,
        normalize_somatic_trees);
 
    // Connect somatic to libraries and save the names of the libraries that
    // were successfully connected.
    first_library_ = num_vertices(pedigree_graph);
    
    add_libraries_to_graph(pedigree_graph, libs);

    num_nodes_ = num_vertices(pedigree_graph);

    // Multiply edge lengths by mutation rates
    update_edge_lengths(pedigree_graph, mu, mu_somatic, mu_library);

    // Remove edges that are non-informative
    simplify_pedigree(pedigree_graph);

    // Prune pedigree
    prune_pedigree(pedigree_graph, inheritance_model_);

    // Apply prefixes to vertex labels to identify germline, somatic, and library nodes
    prefix_vertex_labels(pedigree_graph);

    // Convert vertices in the graph into nodes for peeling operations
    std::vector<size_t> node_ids = ConstructNodes(pedigree_graph);

    family_labels_t family_labels; // (num_families);
    pivots_t pivots;  // (num_families, dummy_index);
    CreateFamiliesInfo(pedigree_graph, &family_labels, &pivots);

    CreatePeelingOps(pedigree_graph, node_ids, family_labels, pivots);
    ConstructPeelingMachine();

    return true;
}

void dng::RelationshipGraph::ConstructPeelingMachine() {
    using namespace dng::peel;
    peeling_functions_.clear();
    peeling_functions_ops_.clear();
    peeling_reverse_functions_.clear();
    peeling_functions_.reserve(peeling_ops_.size());
    peeling_functions_ops_.reserve(peeling_ops_.size());
    peeling_reverse_functions_.reserve(peeling_ops_.size());
    std::vector<std::size_t> lower_written(num_nodes_, -1);
    for(std::size_t i = 0 ; i < peeling_ops_.size(); ++i) {
        peel::Op a = peeling_ops_[i];
        const auto &fam = family_members_[i];
        int b = (int)a;
        auto w = fam[info[b].writes_to];
        bool do_fast = false;
        switch(a) {
        case Op::DOWN:
            // If the lower of the parent has never been written to, we can use the fast version
            do_fast = (lower_written[fam[0]] != -1);
            break;
        case Op::TOCHILD:
            // If the we only have one child, we can use the fast version
            do_fast = (fam.size() == 3);
            break;
        case Op::TOMOTHER:
        case Op::TOFATHER:
        case Op::UP:
            // If the lower of the destination has never been written to, we can use the fast version
            do_fast = (lower_written[w] == -1);
            break;
        default:
            assert(false); // should never get here
            break;
        }
        b = do_fast ? (int)Op::UPFAST + b : b;

        peeling_functions_ops_.push_back(static_cast<peel::Op>(b));
        peeling_functions_.push_back(functions[b]);
        peeling_reverse_functions_.push_back(reverse_functions[b]);

        // If the operation writes to a lower value, make note of it
        if(info[b].writes_lower) {
            lower_written[w] = i;
        }
    }
}

std::vector<std::string> dng::RelationshipGraph::BCFHeaderLines() const {
    using namespace std;
    vector<string> ret = {
        "##META=<ID=OriginalMR,Type=Float,Number=1,Description=\"Mutation rate\">",
        "##META=<ID=FatherMR,Type=Float,Number=1,Description=\"Paternal mutation rate\">",
        "##META=<ID=MotherMR,Type=Float,Number=1,Description=\"Maternal mutation rate\">",
        "##META=<ID=Ploidy,Type=Integer,Number=1,Description=\"Ploidy\">",
        "##META=<ID=Germline,Type=Flag,Number=0,Description=\"Contains germline events\">",
        "##META=<ID=Somatic,Type=Flag,Number=0,Description=\"Contains somatic events\">"
    };

    for(size_t child = 0; child != transitions_.size(); ++child) {
        auto & parents = transitions_[child];
        string line;

        switch(parents.type) {
        case TransitionType::Trio:
            line += "##PEDIGREE=<Child=" + labels_[child];
            line += ",Father=" + labels_[parents.parent1];
            line += ",Mother=" + labels_[parents.parent2];
            line += ",FatherMR=" + utility::to_pretty(parents.length1);
            line += ",MotherMR=" + utility::to_pretty(parents.length2);
            break;
        case TransitionType::Pair:
            line += "##PEDIGREE=<Derived=" + labels_[child];
            line += ",Original=" + labels_[parents.parent1];
            line += ",OriginalMR=" + utility::to_pretty(parents.length1);
            break;
        case TransitionType::Founder:
        default:
            break;    
        }
        line += ",Ploidy=" + utility::to_pretty(ploidies_[child]);
        if(parents.is_germline) {
            line += ",Germline=1";
        }
        if(parents.is_somatic) {
            line += ",Somatic=1";
        }
        if(parents.is_library) {
            line += ",Library=1";
        }        
        ret.push_back(line + ">");
    }
    return ret;
}

namespace {

void prune_pedigree_autosomal(Graph &pedigree_graph);
void prune_pedigree_ylinked(Graph &pedigree_graph);
void prune_pedigree_xlinked(Graph &pedigree_graph);
void prune_pedigree_wlinked(Graph &pedigree_graph);
void prune_pedigree_zlinked(Graph &pedigree_graph);
void prune_pedigree_maternal(Graph &pedigree_graph);
void prune_pedigree_paternal(Graph &pedigree_graph);

void prune_pedigree(Graph &pedigree_graph, InheritanceModel model) {
    switch(model) {
    case InheritanceModel::Autosomal:
        return prune_pedigree_autosomal(pedigree_graph);
    case InheritanceModel::YLinked:
        return prune_pedigree_ylinked(pedigree_graph);
    case InheritanceModel::XLinked:
        return prune_pedigree_xlinked(pedigree_graph);
    case InheritanceModel::WLinked:
        return prune_pedigree_wlinked(pedigree_graph);
    case InheritanceModel::ZLinked:
        return prune_pedigree_zlinked(pedigree_graph);
    case InheritanceModel::Maternal:
        return prune_pedigree_maternal(pedigree_graph);
    case InheritanceModel::Paternal:
        return prune_pedigree_paternal(pedigree_graph);
    default:
        throw std::invalid_argument("Selected inheritance model invalid or not implemented yet.");
    }
}

void prune_pedigree_autosomal(Graph &pedigree_graph) {
    /*do nothing*/;
}

void prune_pedigree_ylinked(Graph &pedigree_graph) {
    auto sexes  = get(boost::vertex_sex, pedigree_graph);
    auto ploidies  = get(boost::vertex_ploidy, pedigree_graph);

    auto is_x = [&](edge_t e) -> bool {
        if(get(boost::edge_type, pedigree_graph, e) == EdgeType::Maternal) {
            return true;
        }
        vertex_t a = source(e, pedigree_graph);
        vertex_t b = target(e, pedigree_graph);
        return (get(boost::vertex_sex, pedigree_graph, std::max(a,b)) != Sex::Male);
    };
    remove_edge_if(is_x, pedigree_graph);

    auto vertex_range = boost::make_iterator_range(vertices(pedigree_graph));
    for (vertex_t v : vertex_range) {
        switch(sexes[v]) {
        case Sex::Female:
            clear_vertex(v,pedigree_graph);
            ploidies[v] = 0;
            break;
        case Sex::Male:
            ploidies[v] = 1;
            break;
        case Sex::Unknown:
        default:
            if(out_degree(v, pedigree_graph) != 0) {
                throw std::invalid_argument("Y-linked inheritance requires every individual to have a known sex.");
            }
        }
    }
}

void prune_pedigree_xlinked(Graph &pedigree_graph) {
    auto is_y = [&](edge_t e) -> bool {
        auto type = get(boost::edge_type, pedigree_graph, e);
        vertex_t a = source(e, pedigree_graph);
        vertex_t b = target(e, pedigree_graph);
        auto sex = get(boost::vertex_sex, pedigree_graph, std::max(a,b));
        auto ploidy = get(boost::vertex_ploidy, pedigree_graph, std::max(a,b));
        if(ploidy == 1) {
            return (sex == Sex::Male);
        }
        return (type == EdgeType::Paternal && sex == Sex::Male);
    };
    remove_edge_if(is_y, pedigree_graph);

    auto sexes  = get(boost::vertex_sex, pedigree_graph);
    auto ploidies  = get(boost::vertex_ploidy, pedigree_graph);
    auto vertex_range = boost::make_iterator_range(vertices(pedigree_graph));
    for (vertex_t v : vertex_range) {
        switch(sexes[v]) {
        case Sex::Female:
            ploidies[v] = (ploidies[v] == 2) ? 2 : 1;
            break;
        case Sex::Male:
            ploidies[v] = 1;
            break;
        case Sex::Unknown:
        default:
            if(out_degree(v, pedigree_graph) != 0) {
                throw std::invalid_argument("X-linked inheritance requires every individual to have a known sex.");
            }
        }
    }
}

void prune_pedigree_wlinked(Graph &pedigree_graph) {
    auto sexes  = get(boost::vertex_sex, pedigree_graph);
    auto ploidies  = get(boost::vertex_ploidy, pedigree_graph);

    auto is_z = [&](edge_t e) -> bool {
        if(get(boost::edge_type, pedigree_graph, e) == EdgeType::Paternal) {
            return true;
        }
        vertex_t a = source(e, pedigree_graph);
        vertex_t b = target(e, pedigree_graph);
        return (get(boost::vertex_sex, pedigree_graph, std::max(a,b)) != Sex::Female);
    };
    remove_edge_if(is_z, pedigree_graph);

    auto vertex_range = boost::make_iterator_range(vertices(pedigree_graph));
    for (vertex_t v : vertex_range) {
        switch(sexes[v]) {
        case Sex::Male:
            clear_vertex(v,pedigree_graph);
            ploidies[v] = 0;
            break;
        case Sex::Female:
            ploidies[v] = 1;
            break;
        case Sex::Unknown:
        default:
            if(out_degree(v, pedigree_graph) != 0) {
                throw std::invalid_argument("W-linked inheritance requires every individual to have a known sex.");
            }
        }
    }
}

void prune_pedigree_zlinked(Graph &pedigree_graph) {
    auto is_w = [&](edge_t e) -> bool {
        auto type = get(boost::edge_type, pedigree_graph, e);
        vertex_t a = source(e, pedigree_graph);
        vertex_t b = target(e, pedigree_graph);
        auto sex = get(boost::vertex_sex, pedigree_graph, std::max(a,b));
        auto ploidy = get(boost::vertex_ploidy, pedigree_graph, std::max(a,b));
        if(ploidy == 1) {
            return (sex == Sex::Female);
        }
        return (type == EdgeType::Maternal && sex == Sex::Female);
    };
    remove_edge_if(is_w, pedigree_graph);

    auto sexes  = get(boost::vertex_sex, pedigree_graph);
    auto ploidies  = get(boost::vertex_ploidy, pedigree_graph);
    auto vertex_range = boost::make_iterator_range(vertices(pedigree_graph));
    for (vertex_t v : vertex_range) {
        switch(sexes[v]) {
        case Sex::Male:
            ploidies[v] = (ploidies[v] == 2) ? 2 : 1;
            break;
        case Sex::Female:
            ploidies[v] = 1;
            break;
        case Sex::Unknown:
        default:
            if(out_degree(v, pedigree_graph) != 0) {
                throw std::invalid_argument("Z-linked inheritance requires every individual to have a known sex.");
            }
        }
    }
}

void prune_pedigree_maternal(Graph &pedigree_graph) {
    auto is_p = [&](edge_t e) -> bool {
        auto edge_type = get(boost::edge_type, pedigree_graph, e);
        return (edge_type == EdgeType::Paternal || edge_type == EdgeType::Spousal);
    };
    remove_edge_if(is_p, pedigree_graph);

    auto ploidies  = get(boost::vertex_ploidy, pedigree_graph);
    auto vertex_range = boost::make_iterator_range(vertices(pedigree_graph));
    for (vertex_t v : vertex_range) {
        ploidies[v] = 1;
    }
}

void prune_pedigree_paternal(Graph &pedigree_graph) {
    auto is_m = [&](edge_t e) -> bool {
        auto edge_type = get(boost::edge_type, pedigree_graph, e);
        return (edge_type == EdgeType::Maternal || edge_type == EdgeType::Spousal);
    };
    remove_edge_if(is_m, pedigree_graph);

    auto ploidies  = get(boost::vertex_ploidy, pedigree_graph);
    auto vertex_range = boost::make_iterator_range(vertices(pedigree_graph));
    for (vertex_t v : vertex_range) {
        ploidies[v] = 1;
    }
}

std::pair<vertex_t,vertex_t> parse_pedigree_table(Graph &pedigree_graph,
        const Pedigree &pedigree, bool normalize_somatic_trees) {

    using namespace std;
    using Sex = dng::Pedigree::Sex;
    using member_t = unsigned int;
    struct child_t {
        member_t id;
        enum struct Type {
            Paternal, Maternal, Clone, Sperm, Egg
        } type;
    };

    // make a copy of names and sexes from pedigree
    vector<string> pedigree_names;
    vector<Sex> pedigree_sexes;
    vector<int> pedigree_ploidies;
    vector<pair<float,float>> pedigree_lengths;

    pedigree_names.reserve(pedigree.NumberOfMembers());
    pedigree_sexes.reserve(pedigree.NumberOfMembers());
    pedigree_ploidies.reserve(pedigree.NumberOfMembers());
    pedigree_lengths.reserve(pedigree.NumberOfMembers());

    auto has_tag = [](const Pedigree::Member &member, const std::string &tag) {
        for(auto &&a : member.tags) {
            if(boost::algorithm::iequals(a, tag)) {
                return true;
            }
        }
        return false;
    };
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

    for(member_t j = 0; j < pedigree.NumberOfMembers(); ++j) {
        auto & member = pedigree.GetMember(j);
        pedigree_names.push_back(member.name);
        pedigree_sexes.push_back(member.sex);
        pedigree_ploidies.push_back(get_ploidy(member));
        pedigree_lengths.push_back(make_pair(member.dad_length.value_or(1.0),
            member.mom_length.value_or(1.0)));
    }

    // identify parents and founders
    vector<member_t> founders;
    vector<vector<child_t>> pedigree_children(pedigree.NumberOfMembers());
    for(member_t j = 0; j < pedigree.NumberOfMembers(); ++j) {
        auto & member = pedigree.GetMember(j);
        if( has_tag(member, "founder") || ( !member.dad && !member.mom ) ) {
            founders.push_back(j);
            continue;
        }
        if( pedigree_ploidies[j] == 0 /* clone */ ) {
            if(member.dad && member.mom) {
                throw std::invalid_argument("Unable to construct graph for pedigree; clone '"
                    + member.name + "' has two parents instead of one.");
            }
            member_t orig_id;
            child_t child{j, child_t::Type::Clone};
            if(member.dad) {
                orig_id = pedigree.LookupMemberPosition(member.dad.get());
            } else {
                orig_id = pedigree.LookupMemberPosition(member.mom.get());
                // move mom's length to the first slot
                pedigree_lengths[j].first = pedigree_lengths[j].second;
            }
            if(orig_id == pedigree.NumberOfMembers()) {
                throw std::invalid_argument("Unable to construct graph for pedigree; the clone parent of '" +
                    member.name + "' is unknown.");
            }
            pedigree_children[orig_id].push_back(child);
            continue;
        }
        if( pedigree_ploidies[j] == 1 /* haploid/gamete */  ) {
            if(member.dad && member.mom) {
                throw std::invalid_argument("Unable to construct graph for pedigree; gamete '"
                    + member.name + "' has two parents instead of one.");
            }
            member_t orig_id;
            child_t child{j};
            if(member.dad) {
                orig_id = pedigree.LookupMemberPosition(member.dad.get());
                child.type = child_t::Type::Sperm;
                if (pedigree_sexes[orig_id] == Sex::Female) {
                    throw std::invalid_argument("Unable to construct graph for pedigree; the father of '" +
                        member.name + "' is female.");
                }
            } else {
                orig_id = pedigree.LookupMemberPosition(member.mom.get());
                child.type = child_t::Type::Egg;
                if (pedigree_sexes[orig_id] == Sex::Male) {
                    throw std::invalid_argument("Unable to construct graph for pedigree; the mother of '" +
                        member.name + "' is male.");
                }
            }
            if(orig_id == pedigree.NumberOfMembers()) {
                throw std::invalid_argument("Unable to construct graph for pedigree; the gamete parent of '" +
                    member.name + "' is unknown.");
            }
            pedigree_children[orig_id].push_back(child);
            continue;
        }

        // Check for valid ids
        if(!member.dad) {
            throw std::invalid_argument("Unable to construct graph for pedigree; the father of '" +
                member.name + "' is unspecified.");
        }
        if(!member.mom) {
            throw std::invalid_argument("Unable to construct graph for pedigree; the mother of '" +
                member.name + "' is unspecified.");
        }
        member_t dad_id = pedigree.LookupMemberPosition(member.dad.get());
        member_t mom_id = pedigree.LookupMemberPosition(member.mom.get());
        if(dad_id == pedigree.NumberOfMembers()) {
            throw std::invalid_argument("Unable to construct graph for pedigree; the father of '" +
                member.name + "' is unknown.");
        }
        if(mom_id == pedigree.NumberOfMembers()) {
            throw std::invalid_argument("Unable to construct graph for pedigree; the mother of '" +
                member.name + "' is unknown.");
        }
        // Check for sanity
        if (pedigree_sexes[dad_id] == Sex::Female) {
            throw std::invalid_argument("Unable to construct graph for pedigree; the father of '" +
                member.name + "' is female.");
        }
        if (pedigree_sexes[mom_id] == Sex::Male ) {
            throw std::invalid_argument("Unable to construct graph for pedigree; the mother of '" +
                member.name + "' is male.");
        }
        // Selfing is not supported
        if (dad_id == mom_id ) {
            throw std::invalid_argument("Unable to construct graph for pedigree; selfing is not supported; "
                "father and mother of '" + member.name + "' are the same.");
        }
        pedigree_children[dad_id].push_back({j,child_t::Type::Paternal});
        pedigree_children[mom_id].push_back({j,child_t::Type::Maternal});
    }
    // Build Graph
    pedigree_graph.clear();
    for(size_t i = 0; i < pedigree_names.size();++i) {
        // Create a germline vertex with an empty label 
        add_vertex({{},VertexType::Germline},pedigree_graph);
    }
    auto labels = get(boost::vertex_label, pedigree_graph);
    auto sexes  = get(boost::vertex_sex, pedigree_graph);
    auto ploidies  = get(boost::vertex_ploidy, pedigree_graph);

    // structure to hold parents on first visit
    vector<std::optional<vertex_t>> touched(pedigree_children.size());
    
    // map id -> vertex_t
    vector<vertex_t> vertices(pedigree_children.size());

    auto process_samples = [&](member_t child_id, vertex_t child_vertex) {
        // Process newick string
        const auto & samples = pedigree.GetMember(child_id).samples;
        std::size_t current_index = num_vertices(pedigree_graph);
        for(auto && sample : samples) {
            if(!parse_newick(sample, child_vertex,
                pedigree_graph, normalize_somatic_trees) ) {
                    throw std::invalid_argument("Unable to parse somatic data for individual '" +
                        pedigree_names[child_id] + "'.");
            }
        }
        // Mark the sex and ploidy of the somatic nodes
        for (vertex_t i = current_index; i < num_vertices(pedigree_graph); ++i) {
            sexes[i] = sexes[child_vertex];
            ploidies[i] = ploidies[child_vertex];
        }        
    };
    
    // Push founders onto a queue
    size_t counter = 0;
    queue<member_t> visited;
    for(auto &&a : founders) {
        labels[counter] = pedigree_names[a];
        sexes[counter] = pedigree_sexes[a];
        ploidies[counter] = pedigree_ploidies[a];
        process_samples(a, counter);

        vertices[a] = counter++;
        visited.push(a);
    }

    // Continue with the rest of the pedigree
    while(!visited.empty()) {
        auto parent_id = visited.front();
        auto parent_vertex = vertices[parent_id];
        visited.pop();
        for(auto child : pedigree_children[parent_id]) {
            vertex_t child_vertex;
            if(child.type == child_t::Type::Clone) {
                child_vertex = counter++;
                vertices[child.id] = child_vertex;
                labels[child_vertex] = pedigree_names[child.id];
                sexes[child_vertex] = sexes[parent_vertex];
                ploidies[child_vertex] = ploidies[parent_vertex];
                add_edge(parent_vertex, child_vertex, {EdgeType::Mitotic, pedigree_lengths[child.id].first}, pedigree_graph);
            } else if(child.type == child_t::Type::Sperm) {
                child_vertex = counter++;
                vertices[child.id] = child_vertex;
                labels[child_vertex] = pedigree_names[child.id];
                sexes[child_vertex] = pedigree_sexes[child.id];
                ploidies[child_vertex] = 1;
                add_edge(parent_vertex, child_vertex, {EdgeType::Paternal, pedigree_lengths[child.id].first}, pedigree_graph);
            } else if(child.type == child_t::Type::Egg) {
                child_vertex = counter++;
                vertices[child.id] = child_vertex;
                labels[child_vertex] = pedigree_names[child.id];
                sexes[child_vertex] = pedigree_sexes[child.id];
                ploidies[child_vertex] = 1;
                add_edge(parent_vertex, child_vertex, {EdgeType::Maternal, pedigree_lengths[child.id].second}, pedigree_graph);
            } else if(!touched[child.id]) {
                // create vertex for child after we have visited both its parents
                touched[child.id] = parent_vertex;
                continue;
            } else {
                // add vertex
                child_vertex = counter++;
                vertices[child.id] = child_vertex;
                labels[child_vertex] = pedigree_names[child.id];
                sexes[child_vertex] = pedigree_sexes[child.id];
                ploidies[child_vertex] = 2;
                // add the meiotic edges
                auto dad_vertex = (child.type != child_t::Type::Maternal) ? parent_vertex : touched[child.id].get();
                auto mom_vertex = (child.type == child_t::Type::Maternal) ? parent_vertex : touched[child.id].get();
                add_edge(mom_vertex, child_vertex, {EdgeType::Maternal, pedigree_lengths[child.id].second}, pedigree_graph);
                add_edge(dad_vertex, child_vertex, {EdgeType::Paternal, pedigree_lengths[child.id].first}, pedigree_graph);

                // Check to see if mom and dad have been seen before
                auto id = edge(dad_vertex, mom_vertex, pedigree_graph);
                if (!id.second) { //Connect dad-mom to make a trio
                    add_edge(dad_vertex, mom_vertex, EdgeType::Spousal, pedigree_graph);
                }
            }
            process_samples(child.id, child_vertex);
            visited.push(child.id);
        }
    }
    // Sanity Check
    assert(counter == pedigree_names.size());
    return {founders.size(),counter};
}

void add_samples_to_graph(Graph &pedigree_graph, const samples_t &known_samples) {
    auto sexes  = get(boost::vertex_sex, pedigree_graph);
    auto ploidies  = get(boost::vertex_ploidy, pedigree_graph);
    auto labels = get(boost::vertex_label, pedigree_graph);
    auto types  = get(boost::vertex_type, pedigree_graph);
    auto library_labels = get(boost::vertex_library_label, pedigree_graph);

    // Create a label-to-id map
    std::map<std::string, vertex_t> soma;
    auto vertex_range = boost::make_iterator_range(vertices(pedigree_graph));
    for(vertex_t v : vertex_range) {
        if(types[v] == VertexType::Somatic && !labels[v].empty()) {
            soma[labels[v]] = v;
        }
    }    

    // Add library nodes to graph
    for (size_t i = 0; i < libs.names.size(); ++i) {
        auto it = soma.find(libs.samples[i]);
        if(it == soma.end()) {
            continue;
        }
        vertex_t u = it->second;
        std::string name = labels[u];
        if(libs.names[i] != labels[u]) {
            name += DNG_LABEL_SEPARATOR;
            name += libs.names[i];
        }
        vertex_t v = add_vertex({name, VertexType::Library}, pedigree_graph);
        add_edge(u, v, {EdgeType::Library, 1.0f}, pedigree_graph);
        sexes[v] = sexes[u];
        ploidies[v] = ploidies[u];
        library_labels[v] = libs.names[i];
    }
}

void update_edge_lengths(Graph &pedigree_graph,
        double mu_meiotic, double mu_somatic, double mu_library) {
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    auto edge_types = get(boost::edge_type, pedigree_graph);
    auto lengths = get(boost::edge_length, pedigree_graph);

    auto range = boost::make_iterator_range(edges(pedigree_graph));
    for (edge_t e : range) {
        switch(edge_types[e]) {
        case EdgeType::Maternal:
        case EdgeType::Paternal:
            lengths[e] *= mu_meiotic;
            break;
        case EdgeType::Mitotic:
            lengths[e] *= mu_somatic;
            break;
        case EdgeType::Library:
            lengths[e] *= mu_library;
            break;
        default:
            break;
        }
    }
}

void prefix_vertex_labels(dng::detail::graph::Graph &pedigree_graph) {
    using namespace dng::detail::graph;

    auto labels = get(boost::vertex_label, pedigree_graph);
    auto types = get(boost::vertex_type, pedigree_graph);

    auto vertex_range = boost::make_iterator_range(vertices(pedigree_graph));
    for(vertex_t v : vertex_range) {
        const char *ch = "";
        switch(types[v]) {
        case VertexType::Germline:
            ch = (DNG_LABEL_PREFIX_GERMLINE DNG_LABEL_SEPARATOR);
            break;
        case VertexType::Somatic:
            ch = (DNG_LABEL_PREFIX_SOMATIC DNG_LABEL_SEPARATOR);
            break;
        case VertexType::Library:
            ch = (DNG_LABEL_PREFIX_LIBRARY DNG_LABEL_SEPARATOR);
            break;
        default:
            break;
        }
        labels[v].insert(0,ch);
    }
}


void simplify_pedigree(Graph &pedigree_graph) {
    using boost::make_iterator_range;

    auto edge_types = get(boost::edge_type, pedigree_graph);
    auto lengths = get(boost::edge_length, pedigree_graph);
    auto types = get(boost::vertex_type, pedigree_graph);

    auto vertex_range = boost::adaptors::reverse(make_iterator_range(vertices(pedigree_graph)));

    for (vertex_t v : vertex_range) {
        if(types[v] == VertexType::Library) {
            continue;
        }
        // identify children
        auto edge_range = make_iterator_range(out_edges(v, pedigree_graph));
        // NOTE: This assumes that children is a reasonable number
        // and will kill performance if it is not.
        std::vector<vertex_t> meiotic_children;
        for(edge_t && e : edge_range) {
            vertex_t u = target(e, pedigree_graph);
            if( u > v &&
                (edge_types[e] == EdgeType::Maternal || 
                 edge_types[e] == EdgeType::Paternal)) {
                meiotic_children.push_back(u);
            }
        }

        // Remove spousal edges with no children
        remove_out_edge_if(v, [&](edge_t e) -> bool {
            if(edge_types[e] != EdgeType::Spousal) {
                return false;
            }
            vertex_t u = target(e, pedigree_graph);
            auto spousal_edge_range = make_iterator_range(out_edges(u, pedigree_graph));
            for(edge_t && e : spousal_edge_range) {
                vertex_t o = target(e, pedigree_graph);
                auto it = boost::range::find(meiotic_children, o);
                if(it != meiotic_children.end()) {
                    return false;
                }
            }
            return true;
        }, pedigree_graph);

        size_t children = 0, ancestors = 0, spouses = 0;
        edge_range = make_iterator_range(out_edges(v, pedigree_graph));
        for (auto && e : edge_range) {
            if (edge_types[e] == EdgeType::Spousal) {
                spouses += 1;
            } else if (target(e, pedigree_graph) > v) {
                children += 1;
            } else {
                ancestors += 1;
            }
        }
        if (children == 0) {
            // this node has no descendants
            clear_vertex(v, pedigree_graph);
        } else if (children >= 2 || spouses != 0) {
            /*noop*/;
        } else if (ancestors > 0) {
            assert(ancestors < 3); 
            edge_t edge_trio[3];
            vertex_t vertex_trio[3];//
            int child_index = ancestors;
            auto it = edge_range.begin();
            for (size_t j = 0; j <= ancestors; ++j) {
                edge_trio[j] = *it++;
                vertex_trio[j] = target(edge_trio[j], pedigree_graph);
            }
            for (size_t p = 0; p < child_index; ++p) {
                if (vertex_trio[child_index] < vertex_trio[p]) {
                    boost::swap(vertex_trio[child_index], vertex_trio[p]);
                    boost::swap(edge_trio[child_index], edge_trio[p]);
                }
                add_edge(vertex_trio[p], vertex_trio[child_index],
                         {edge_types[edge_trio[p]], lengths[edge_trio[p]] + lengths[edge_trio[child_index]]},
                         pedigree_graph);
            }
            clear_vertex(v, pedigree_graph);
        }
    }
}

} // namespace 

std::vector<size_t> dng::RelationshipGraph::ConstructNodes(const Graph &pedigree_graph) {
    // node_ids[vertex] will convert vertex_id to node position.
    std::vector<size_t> node_ids(num_vertices(pedigree_graph), -1);

    labels_.clear();
    labels_.reserve(128);
    ploidies_.clear();
    ploidies_.reserve(128);

    auto labels = get(boost::vertex_label, pedigree_graph);
    auto ploidies = get(boost::vertex_ploidy, pedigree_graph);
    auto library_labels = get(boost::vertex_library_label, pedigree_graph);

    auto vertex_range = boost::make_iterator_range(vertices(pedigree_graph));
    for(vertex_t v : vertex_range) {
        // Skip vertices with no edges
        if (out_degree(v, pedigree_graph) == 0) {
            continue;
        }
        if(!library_labels[v].empty()) {
            library_names_.push_back(library_labels[v]);
        }
        const auto vid = labels_.size();
        node_ids[v] = vid;

        if (labels[v].empty() || labels[v].back() == DNG_LABEL_SEPARATOR_CHAR) {
            labels_.push_back(labels[v] + "unnamed_node_" + utility::to_pretty(vid));
        } else {
            labels_.push_back(labels[v]);
        }
        ploidies_.push_back(ploidies[v]);
    }
    num_nodes_ = labels_.size();

    auto update_position = [&node_ids](size_t pos, size_t last) -> size_t {
        for (; pos < node_ids.size() && node_ids[pos] == -1; ++pos)
            /*noop*/;
        return (pos < node_ids.size()) ? node_ids[pos] : last;
    };

    first_founder_ = update_position(first_founder_, num_nodes_);
    first_nonfounder_ = update_position(first_nonfounder_, num_nodes_);
    first_somatic_ = update_position(first_somatic_, num_nodes_);
    first_library_ = update_position(first_library_, num_nodes_);

    return node_ids;
}

void dng::RelationshipGraph::CreateFamiliesInfo(Graph &pedigree_graph,
        family_labels_t *family_labels, pivots_t *pivots)
{
    assert(family_labels != nullptr);
    assert(pivots != nullptr);

    auto groups = get(boost::vertex_group, pedigree_graph);
    auto families = get(boost::edge_family, pedigree_graph);

    // Calculate the connected components.  This defines independent sections
    // of the graph.
    std::size_t num_groups = connected_components(pedigree_graph, groups);

    // Calculate the biconnected components and articulation points.
    // This defines "nuclear" families and pivot individuals.
    // Nodes which have no edges will not be part of any family.
    std::vector<vertex_t> articulation_vertices;
    std::size_t num_families = biconnected_components(pedigree_graph, families,
                                   back_inserter(articulation_vertices)).first;

    family_labels->assign(num_families, {});
    pivots->assign(num_families, {});

    // Determine which edges belong to which nuclear families.
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
        (*family_labels)[families[*ei]].push_back(*ei);
    }

    // Determine the last family in each group.  All singleton groups will have
    // a value of -1 since they have no family assignment.
    using root_families_t = std::vector<boost::optional<std::size_t>>;
    root_families_t root_families(num_groups);
    for (std::size_t f = 0; f < family_labels->size(); ++f) {
        // last one wins
        auto first_edge = (*family_labels)[f][0];
        auto src_vertex = source(first_edge, pedigree_graph);
        root_families[groups[src_vertex]] = f;
    }

    // Identify the pivot for each family.
    // The pivot will be the last art. point that has an edge in
    // the group.  The pivot of the last group doesn't matter.
    for (auto a : articulation_vertices) {
        boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = out_edges(a, pedigree_graph); ei != ei_end; ++ei) {
            // Just overwrite existing value so that the last one wins.
            (*pivots)[families[*ei]] = a;
        }
    }

    // Clear pivots of roots
    for (auto f : root_families) {
        if(f) {
            (*pivots)[f.get()] = boost::none;
        }
    }
}

void dng::RelationshipGraph::CreatePeelingOps(
        const Graph &pedigree_graph, const std::vector<size_t> &node_ids,
        family_labels_t &family_labels, const pivots_t &pivots) {

    auto edge_types = get(boost::edge_type, pedigree_graph);
    auto lengths = get(boost::edge_length, pedigree_graph);
    auto vertex_types = get(boost::vertex_type, pedigree_graph);

    using VertexType = detail::graph::VertexType;

    ClearFamilyInfo();

    transitions_.resize(num_nodes_);

    constexpr size_t null_id = static_cast<size_t>(-1);

    // Setup founder Transitions
    for (std::size_t i = first_founder_; i < first_nonfounder_; ++i) {
        transitions_[i] = {TransitionType::Founder, null_id, null_id, 0.0, 0.0,false,false,false};
    }

    // Detect Family Structure and pivot positions
    for (std::size_t k = 0; k < family_labels.size(); ++k) {
        auto &family_edges = family_labels[k];

        // Sort edges based on type and target
        boost::sort(family_edges, [&](edge_t x, edge_t y) -> bool {
            if(edge_types(x) < edge_types(y)) {
                return true;
            }
            if(edge_types(x) == edge_types(y)) {
                return target(x, pedigree_graph) < target(y, pedigree_graph);   
            }
            return false;
        });

        // Find the range of the parent types
        auto it = boost::find_if(family_edges, [&](edge_t x) -> bool {
            return (edge_types(x) != EdgeType::Spousal);
        });
        size_t num_spousal_edges = std::distance(family_edges.begin(), it);

        // Check to see what type of graph we have
        if (num_spousal_edges == 0) {
            // If we do not have a parent-child single branch,
            // we can't construct the pedigree.
            if (family_edges.size() != 1) {
                throw std::invalid_argument("Unable to construct peeler for pedigree;  "
                        "it does not have a parent-child single branch");
            }
            // Create a mitotic peeling operation.
            auto child_index = target(*it, pedigree_graph);
            auto parent_index = source(*it, pedigree_graph);
            size_t parent = node_ids[parent_index];
            size_t child = node_ids[child_index];
            
            transitions_[child] = {TransitionType::Pair, parent, null_id, lengths[*it], 0.0};

            auto child_type = vertex_types(child_index);
            auto parent_type = vertex_types(parent_index);

            bool temp_types[3] = {false,false,false};
            for(int x = static_cast<int>(parent_type); x <= static_cast<int>(child_type); ++x) {
                temp_types[x] = true;
            }
            transitions_[child].is_germline = temp_types[0];
            transitions_[child].is_somatic  = temp_types[1];
            transitions_[child].is_library  = temp_types[2];

            family_members_.push_back({parent, child});

            if(pivots[k]) {
                if(node_ids[pivots[k].get()] == child) {
                    peeling_ops_.push_back(peel::Op::DOWN);
                } else {
                    peeling_ops_.push_back(peel::Op::UP);
                }                
            } else {
                peeling_ops_.push_back(peel::Op::UP);
                roots_.push_back(parent);
            }
        } else if(num_spousal_edges == 1) {
            // If this family contains no children, skip it
            if (it == family_edges.end()) {
                continue;
            }
            // We have a nuclear family with 1 or more children

            auto dad_index = source(family_edges.front(), pedigree_graph);
            auto mom_index = target(family_edges.front(), pedigree_graph);
            size_t dad = node_ids[dad_index];
            size_t mom = node_ids[mom_index];

            family_members_.push_back({dad, mom});
            auto &family_members = family_members_.back();

            for (;it != family_edges.end();++it) {
                auto child_index = target(*it, pedigree_graph);

                size_t child = node_ids[child_index];
                if(edge_types(*it) == EdgeType::Maternal) {
                    transitions_[child] = {TransitionType::Trio, dad, mom, 0, lengths[*it]};
                    
                    auto child_type = vertex_types(child_index);
                    auto parent_type = vertex_types(mom_index);

                    bool temp_types[3] = {false,false,false};
                    for(int x = static_cast<int>(parent_type); x <= static_cast<int>(child_type); ++x) {
                        temp_types[x] = true;
                    }
                    transitions_[child].is_germline = temp_types[0];
                    transitions_[child].is_somatic  = temp_types[1];
                    transitions_[child].is_library  = temp_types[2];

                    family_members.push_back(child);
                } else {
                    assert(edge_types(*it) == EdgeType::Paternal);
                    assert(transitions_[child].type == TransitionType::Trio);
                    transitions_[child].length1 = lengths[*it];
                }
            }
            if (!pivots[k]) {
                // A family without a pivot is a root family
                peeling_ops_.push_back(peel::Op::TOFATHER);
                roots_.push_back(family_members[0]);

            } else {
                auto pivot_pos = boost::range::find(family_members, node_ids[pivots[k].get()]);
                size_t p = distance(family_members.begin(), pivot_pos);

                if (p == 0) {
                    peeling_ops_.push_back(peel::Op::TOFATHER);
                } else if (p == 1) {
                    peeling_ops_.push_back(peel::Op::TOMOTHER);
                } else if (p == 2) {
                    peeling_ops_.push_back(peel::Op::TOCHILD);
                } else {
                    peeling_ops_.push_back(peel::Op::TOCHILD);
                    boost::swap(family_members[p], family_members[2]);
                }
            }

        } else {
            throw std::invalid_argument("Unable to construct peeler for pedigree; Not a zero-loop pedigree");
        }
    }
}

void dng::RelationshipGraph::ClearFamilyInfo(){
    roots_.clear();
    roots_.reserve(16);
    family_members_.clear();
    family_members_.reserve(128);
    peeling_ops_.clear();
    peeling_ops_.reserve(128);
    transitions_.clear();
    transitions_.reserve(128);
}
