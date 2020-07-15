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
#include <mutk/mutation.hpp>

#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/fill.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/range/numeric.hpp>

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

struct workspace_t {
    std::vector<mutk::Tensor<1>> stack;
    std::array<mutk::tensor_index_t,3> widths;
    float scale;
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

    float PeelForward(workspace_t *work) const;

    void PrintGraph(std::ostream &os) const;

    samples_t SampleNames() const;

    workspace_t CreateWorkspace() const;

    const std::vector<potential_t>& potentials() const {
        return potentials_;
    }

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

    // Elimination Rank
    std::vector<std::size_t> elimination_rank_;

    // Peeling Operations
    using peeling_t = std::unique_ptr<PeelingVertex>;
    std::vector<peeling_t> peelers_;

    // location of root data in peelers
    std::vector<size_t> roots_;

private:

};

class PeelingVertex {
public:
    using workspace_t = mutk::RelationshipGraph::workspace_t;
    using tree_t = detail::junction_tree::Graph;
    using vertex_t = detail::junction_tree::vertex_t;
    using label_t = mutk::detail::pedigree_graph::vertex_t;

    virtual ~PeelingVertex() = default;

    virtual void Forward(PeelingVertex::workspace_t *work) const = 0;

    virtual void AddLocal(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) = 0;
    virtual void AddOutput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) = 0;
    virtual void AddInput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) = 0;
};

template<std::size_t N>
class GeneralPeelingVertex : public PeelingVertex {
public:
    using tensor_t = mutk::Tensor<N>;

    struct data_t {
        std::size_t index;
        typename tensor_t::Dimensions data_ploidy;
        typename tensor_t::Dimensions broadcast_ploidy;
    };

    GeneralPeelingVertex(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index, const std::vector<int> &buffer_ploidy);

    virtual void AddLocal(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) override;
    virtual void AddOutput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) override;
    virtual void AddInput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) override;

    virtual ~GeneralPeelingVertex() = default;

    virtual void Forward(PeelingVertex::workspace_t *work) const override;

protected:
    data_t MakeMetadata(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) const;

    std::array<PeelingVertex::label_t,N> labels_;

    data_t buffer_;

    data_t local_data_;
    data_t output_data_;
    std::vector<data_t> input_data_;
};

namespace detail {

template<typename A>
auto calc_dims(const PeelingVertex::workspace_t &work, A a) {
    typename mutk::Tensor<a.size()>::Dimensions ret;
    std::transform(a.begin(), a.end(), ret.begin(), [&](auto v) {
        return work.widths[v];
    });
    return ret;
}

}

template<std::size_t N>
GeneralPeelingVertex<N>::GeneralPeelingVertex(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index, const std::vector<int> &buffer_ploidy) {
    boost::range::copy(labels, labels_.begin());
    buffer_.index = buffer_index;
    boost::range::copy(buffer_ploidy, buffer_.data_ploidy.begin());
    boost::range::fill(buffer_.broadcast_ploidy, 0);
}

template<std::size_t N>
typename GeneralPeelingVertex<N>::data_t
GeneralPeelingVertex<N>::MakeMetadata(const std::vector<PeelingVertex::label_t> &labels,
    std::size_t buffer_index) const {
    data_t ret;
    ret.index = buffer_index;

    // Go through the labels for this vertex and see if they are in the Metadata labels
    for(size_t i = 0; i < labels_.size(); ++i) {
        auto it = boost::range::find(labels, labels_[i]);
        if(it != labels.end()) {
            // our metadata contains the label, we will not need to broadcast it
            ret.data_ploidy[i] = buffer_.data_ploidy[i];
            ret.broadcast_ploidy[i] = 0;
        } else {
            // we will need to broadcast this axis
            ret.data_ploidy[i] = 0;
            ret.broadcast_ploidy[i] = buffer_.data_ploidy[i];
        }
    }
    return ret;
}

template<std::size_t N>
void GeneralPeelingVertex<N>::AddLocal(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) {
    local_data_ = MakeMetadata(labels, buffer_index);
}

template<std::size_t N>
void GeneralPeelingVertex<N>::AddOutput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) {
    output_data_ = MakeMetadata(labels, buffer_index);
}

template<std::size_t N>
void GeneralPeelingVertex<N>::AddInput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) {
    input_data_.push_back(MakeMetadata(labels, buffer_index));
}

template<std::size_t N>
void GeneralPeelingVertex<N>::Forward(PeelingVertex::workspace_t *work) const {
    // copy local data to local buffer
    auto dims = detail::calc_dims(*work, local_data_.data_ploidy);
    if(work->stack[local_data_.index].size() == 0) {
        // unit potentials are empty
        work->stack[buffer_.index].resize(dims.TotalSize());
        work->stack[buffer_.index].setConstant(1.0f);
    } else {
        work->stack[buffer_.index] = work->stack[local_data_.index];
    }

    for(int i=0;i<input_data_.size(); ++i) {
        auto input_dims = detail::calc_dims(*work, input_data_[i].data_ploidy);
        auto broadcast_dims = detail::calc_dims(*work, input_data_[i].broadcast_ploidy);
        work->stack[buffer_.index].reshape(dims) *= 
            work->stack[input_data_[i].index].reshape(input_dims).broadcast(broadcast_dims);
    }
    auto msg_dims = detail::calc_dims(*work, output_data_.data_ploidy);
    std::vector<typename tensor_t::Index> sum_dims;

    for(typename tensor_t::Index i=0; i < output_data_.data_ploidy.size(); ++i) {
        if(output_data_.data_ploidy[i] == 0) {
            sum_dims.push_back(i);
        }
    }

    work->stack[output_data_.index].resize(msg_dims.TotalSize());
    switch(sum_dims.size()) {
    case 0:
        work->stack[output_data_.index] = work->stack[buffer_.index];
        break;
    case 1:
        if constexpr( N >= 1 ) {
            work->stack[output_data_.index].reshape(msg_dims) =
                work->stack[buffer_.index].reshape(dims)
                    .sum(Eigen::make_index_list(sum_dims[0])).reshape(msg_dims);
            break;
        }
    case 2:
        if constexpr( N >= 2 ) {
            work->stack[output_data_.index].reshape(msg_dims) =
                work->stack[buffer_.index].reshape(dims)
                    .sum(Eigen::make_index_list(sum_dims[0],sum_dims[1])).reshape(msg_dims);
            break;            
        }
    case 3:
        if constexpr( N >= 3 ) {
            work->stack[output_data_.index].reshape(msg_dims) =
                work->stack[buffer_.index].reshape(dims)
                    .sum(Eigen::make_index_list(sum_dims[0],sum_dims[1],sum_dims[2])).reshape(msg_dims);            
            break;
        }
    case 4:
        if constexpr( N >= 4 ) {
            work->stack[output_data_.index].reshape(msg_dims) =
                work->stack[buffer_.index].reshape(dims)
                    .sum(Eigen::make_index_list(sum_dims[0],sum_dims[1],sum_dims[2],sum_dims[3])).reshape(msg_dims);            
            break;
        }        
    default:
        assert(false);
    };
}

}; // namespace mutk

#endif // MUTK_RELATIONSHIP_GRAPH_HPP
