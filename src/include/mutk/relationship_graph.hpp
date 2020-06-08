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

#ifndef MUTK_RELATIONSHIP_GRAPH_HPP
#define MUTK_RELATIONSHIP_GRAPH_HPP

#include <mutk/detail/graph.hpp>
#include <mutk/pedigree.hpp>
#include <mutk/memory.hpp>
#include <mutk/peeling.hpp>

#include <boost/container/static_vector.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace mutk {

enum struct InheritanceModel {
    Autosomal,    // default option
    Maternal,     // transmitted by mother to child
    Paternal,     // transmitter by father to child
    XLinked,      // females have 2 copies, males have 1; males transmit to daughters, not to sons
    YLinked,      // males have 1 copy, only transmits it to sons
    WLinked,      // females have 1 copy, only transmitted to daughters
    ZLinked,      // males have 2 copies, females have 1; females transmit to sons, not to daughters
    Invalid
};

namespace detail {
extern const std::map<std::string, InheritanceModel> CHR_MODEL_MAP;

enum struct Potential {
    LikelihoodDiploid,   // P(Data|G)
    LikelihoodHaploid,   // P(Data|H)
    FounderDiploid,      // P(G)
    FounderHaploid,      // P(H)
    CloneDiploid,        // P(G1|G2)
    CloneHaploid,        // P(H1|H2)
    GameteDiploid,       // P(H1|G2)
    ChildDiploidDiploid, // P(G1|G2,G3)
    ChildHaploidDiploid, // P(G1|H2,G3)
    ChildDiploidHaploid, // P(G1|G2,H3)
    ChildHaploidHaploid, // P(G1|H1,H2)
    ChildSelfingDiploid, // P(G1|G2,G2)
    ChildSelfingHaploid  // P(G1|H2,H2)
};

struct potential_t {
    Potential type;
    int child;
    using data_t = std::pair<int,float>;
    boost::container::static_vector<data_t, 2> parents;

    potential_t(Potential type_arg, int child_arg) : type{type_arg}, child{child_arg} {}
    potential_t(Potential type_arg, int child_arg, int par1, float dist1) :
        type{type_arg}, child{child_arg}, parents{{par1,dist1}} {}
    potential_t(Potential type_arg, int child_arg, int par1, float dist1, int par2, float dist2) :
        type{type_arg}, child{child_arg}, parents{{par1,dist1},{par2,dist2}} {}
};

struct workspace_t {
    std::vector<mutk::matrix_t> stack;
};

} //namespace detail

class RelationshipGraph {
public:
    using samples_t = std::vector<const char*>;
    using workspace_t = detail::workspace_t;
    using Potential = detail::Potential;
    using potential_t = detail::potential_t;

    void ConstructGraph(const Pedigree& pedigree, const samples_t& known_samples,
            InheritanceModel inheritance_model,
            double mu, double mu_somatic, bool normalize_somatic_trees);

    void ConstructMachine();


    void PrintGraph(std::ostream &os) const;

    samples_t SampleNames() const;

    // double PeelForwards(peel::workspace_t &work,
    //                     const TransitionMatrixVector &mat) const {
    //     if(work.dirty_lower) {
    //         work.CleanupFast();
    //     }

    //     // Peel pedigree one family at a time
    //     for(std::size_t i = 0; i < peeling_functions_.size(); ++i) {
    //         (*peeling_functions_[i])(work, family_members_[i], mat);
    //     }

    //     // Sum over roots
    //     double ret = 0.0;
    //     for(auto r : roots_) {
    //         ret += log((work.lower[r] * work.upper[r]).sum());
    //     }
        
    //     return ret;
    // }

    // double PeelBackwards(peel::workspace_t &work,
    //                      const TransitionMatrixVector &mat) const {
    //     double ret = 0.0;
    //     // Divide by the likelihood
    //     for(auto r : roots_) {
    //         double sum = (work.lower[r] * work.upper[r]).sum();
    //         ret += log(sum);
    //         work.lower[r] /= sum;
    //     }

    //     for(std::size_t i = peeling_reverse_functions_.size(); i > 0; --i) {
    //         (*peeling_reverse_functions_[i - 1])(work, family_members_[i - 1], mat);
    //     }
    //     work.dirty_lower = true;
    //     return ret;
    // }

    workspace_t CreateWorkspace() const;

protected:
    InheritanceModel inheritance_model_{InheritanceModel::Autosomal};

    detail::pedigree_graph::Graph graph_;
    detail::junction_tree::Graph junction_tree_;

    using vertex_range_t = std::pair<detail::pedigree_graph::vertex_t,detail::pedigree_graph::vertex_t>;

    vertex_range_t founders_;
    vertex_range_t descendants_;
    vertex_range_t samples_;

    size_t stack_size_;

    // The probability potentials of the relationship graph
    std::vector<potential_t> potentials_;

private:


};

// template<typename A, typename M>
// inline
// RelationshipGraph create_relationship_graph(const A& arg, M* mpileup) {
//     assert(mpileup != nullptr);
//     // Parse Pedigree from File
//     Pedigree ped = io::parse_ped(arg.ped);
//     // Construct peeling algorithm from parameters and pedigree information
//     RelationshipGraph relationship_graph;
//     if (!relationship_graph.Construct(ped, mpileup->libraries(), inheritance_model(arg.model),
//                                       arg.mu, arg.mu_somatic, arg.mu_library,
//                                       arg.normalize_somatic_trees)) {
//         throw std::runtime_error("Unable to construct peeler for pedigree; "
//                                  "possible non-zero-loop relationship_graph.");
//     }
//     // Select libraries in the input that are used in the pedigree
//     mpileup->SelectLibraries(relationship_graph.library_names());

//     return relationship_graph;
// }

}; // namespace mutk

#endif // MUTK_RELATIONSHIP_GRAPH_HPP
