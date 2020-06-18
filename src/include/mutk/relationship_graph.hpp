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
    std::vector<mutk::Tensor<1>> stack;
    std::array<std::size_t,3> widths;
};

} //namespace detail


class PeelingVertex;

class RelationshipGraph {
public:
    using samples_t = std::vector<const char*>;
    using workspace_t = detail::workspace_t;
    using Potential = detail::Potential;
    using potential_t = detail::potential_t;

    using graph_t = detail::pedigree_graph::Graph;
    using junction_tree_t = detail::junction_tree::Graph;

    void ConstructGraph(const Pedigree& pedigree, const samples_t& known_samples,
            InheritanceModel inheritance_model,
            double mu, double mu_somatic, bool normalize_somatic_trees);

    void ConstructPeeler();


    void PrintGraph(std::ostream &os) const;

    samples_t SampleNames() const;

    workspace_t CreateWorkspace() const;

protected:
    InheritanceModel inheritance_model_{InheritanceModel::Autosomal};

    graph_t graph_;
    junction_tree_t junction_tree_;

    using vertex_range_t = std::pair<detail::pedigree_graph::vertex_t,detail::pedigree_graph::vertex_t>;

    vertex_range_t founders_;
    vertex_range_t descendants_;
    vertex_range_t samples_;

    size_t stack_size_;

    // The probability potentials of the relationship graph
    std::vector<potential_t> potentials_;

    using peeling_t = std::unique_ptr<PeelingVertex>;
    std::vector<peeling_t> peeling_vertexes_;

private:

};

class PeelingVertex {
public:
    using workspace_t = mutk::RelationshipGraph::workspace_t;
    using tree_t = detail::junction_tree::Graph;
    using vertex_t = detail::junction_tree::vertex_t;

    virtual ~PeelingVertex() = default;
    virtual void Forward(RelationshipGraph::workspace_t *work) const = 0;
};

template<std::size_t N>
class GeneralPeelingVertex {
public:
    using tensor_t = mutk::Tensor<N>;

    virtual ~GeneralPeelingVertex() = default;

    struct data_t {
        std::size_t index;
        typename tensor_t::Dimensions data_ploidy;
        typename tensor_t::Dimensions broadcast_ploidy;
    };

    void Forward(PeelingVertex::workspace_t *work);

    data_t buffer_;

    data_t local_data_;
    data_t output_data_;
    std::vector<data_t> input_data_;
};

namespace detail {
template<std::size_t N>
using dims_t = typename mutk::Tensor<N>::Dimensions;
}

template<std::size_t N>
void GeneralPeelingVertex<N>::Forward(PeelingVertex::workspace_t *work) {
    // copy local data to local buffer
    work->stack[buffer_.index] = work->stack[local_data_.index];

    auto dims = calc_dims(*work, local_data_.data_ploidy);
    for(int i=0;i<input_data_.size(); ++i) {
        auto input_dims = calc_dims(*work, input_data_[i].data_ploidy);
        auto broadcast_dims = calc_dims(*work, input_data_[i].broadcast_ploidy);
        work->stack[buffer_.index].reshape(dims) *= 
            work->stack[input_data_[i].index].reshape(input_dims).broadcast(broadcast_dims);
    }
    auto msg_dims = calc_dims(*work, output_data_.data_ploidy);
    std::vector<typename tensor_t::Index> sum_dims;

    for(typename tensor_t::Index i=0; i < msg_dims.size(); ++i) {
        if(msg_dims[i] == 1) {
            sum_dims.push_back(i);
        }
    }
    switch(sum_dims.size()) {
    case 0: {
        work->stack[output_data_.index] = work->stack[buffer_.index];
        break;
    }
    case 1: {
        work->stack[output_data_.index].reshape(msg_dims) =
            work->stack[buffer_.index].reshape(dims)
                .sum(detail::dims_t<1>{sum_dims[0]}).reshape(msg_dims);
        break;
    }
    case 2: {
        work->stack[output_data_.index].reshape(msg_dims) =
            work->stack[buffer_.index].reshape(dims)
                .sum(detail::dims_t<2>{sum_dims[0],sum_dims[1]}).reshape(msg_dims);
        break;
    }
    case 3: {
        work->stack[output_data_.index].reshape(msg_dims) =
            work->stack[buffer_.index].reshape(dims)
                .sum(detail::dims_t<3>{sum_dims[0],sum_dims[1],sum_dims[2]}).reshape(msg_dims);
        break;
    }    
    default:
        assert(false);
    };
}



}; // namespace mutk

#endif // MUTK_RELATIONSHIP_GRAPH_HPP
