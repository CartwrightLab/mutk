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
    std::vector<mutk::tensor_t> stack;
    float_t scale;
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

    float PeelForward(workspace_t &work) const;

    void PrintGraph(std::ostream &os) const;

    samples_t SampleNames() const;

    workspace_t CreateWorkspace() const;

    const std::vector<potential_t>& potentials() const {
        return potentials_;
    }

    const graph_t & graph() const {
        return graph_;
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
    std::vector<mutk::shape_t> ploidy_shapes_;
    

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

    virtual void Forward(PeelingVertex::workspace_t &work) const = 0;

    virtual void AddLocal(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) = 0;
    virtual void AddOutput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) = 0;
    virtual void AddInput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) = 0;
};

class GeneralPeelingVertex : public PeelingVertex {
public:

    struct data_t {
        std::size_t index;
        shape_t shape;
    };

    GeneralPeelingVertex(const std::vector<PeelingVertex::label_t> &labels) :
        labels_{labels} { }

    virtual void AddLocal(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) override;
    virtual void AddOutput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) override;
    virtual void AddInput(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) override;

    virtual ~GeneralPeelingVertex() = default;

    virtual void Forward(PeelingVertex::workspace_t &work) const override;

protected:
    data_t MakeMetadata(const std::vector<PeelingVertex::label_t> &labels,
        std::size_t buffer_index) const;

    std::vector<PeelingVertex::label_t> labels_;

    mutable tensor_t temporary_;

    data_t local_data_;
    data_t output_data_;
    std::vector<data_t> input_data_;
};

}

}; // namespace mutk

#endif // MUTK_RELATIONSHIP_GRAPH_HPP
