/*
# Copyright (c) 2014-2020 Reed A. Cartwright <reed@cartwright.ht>
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

#ifndef MUTK_POTENTIAL_HPP
#define MUTK_POTENTIAL_HPP

#include "message.hpp"
#include "mutation.hpp"

#include <boost/container/flat_set.hpp>
#include <vector>

namespace mutk {

enum struct PotentialType {
    Unit,                // All ones
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

class Potential {
 public:
    using labels_t = boost::container::flat_set<message_label_t>;

    explicit Potential(labels_t labels) : labels_{std::move(labels)} {}

    template<class It>
    Potential(It first, It last) : labels_{first, last} {}

    explicit Potential(const std::vector<message_label_t> &seq) : Potential(seq.begin(), seq.end()) {}

    explicit Potential(std::initializer_list<message_label_t> init) : labels_{init} {}

    virtual ~Potential() = default;

    enum struct any_t  : int {};
    enum struct mean_t : int {};
    enum struct some_t : int {};

    virtual message_t Create(message_size_t n, any_t) = 0;
    virtual message_t Create(message_size_t n, some_t) = 0;
    virtual message_t Create(message_size_t n, mean_t) = 0;

    static constexpr any_t  ANY{0};
    static constexpr mean_t MEAN{1};
    static constexpr some_t ZERO{0};
    static constexpr some_t ONE{1};
    static constexpr some_t TWO{2};

   inline
   message_t::shape_type Shape(int n) const {
       message_t::shape_type ret(labels_.size());
       std::transform(labels_.begin(), labels_.end(), ret.begin(),
           [n](auto v) { return message_axis_size(n,v); });
       return ret;
   }

   inline
   std::vector<int> Ploidies() const {
      std::vector<int> ret;
      for( auto && a : labels_) {
         ret.push_back((int)message_axis_ploidy(a));
      }
      return ret;
   }

 protected:
    labels_t labels_;
};

class UnitPotential : public Potential {
 public:
    using Potential::Potential;

    virtual message_t Create(message_size_t n, any_t) override;
    virtual message_t Create(message_size_t n, some_t) override;
    virtual message_t Create(message_size_t n, mean_t) override;
};

class MutationPotential : public Potential {
 public:
    template<typename... Args>
    MutationPotential(const MutationModel & model, Args... args) : Potential(std::forward<Args>(args)...),
        model_(model),
        probability_builder_(Ploidies()), expectation_builder_(Ploidies()),
        zero_count_builder_(Ploidies()), one_count_builder_(Ploidies())
    { }

    virtual message_t Create(message_size_t n, any_t) override;
    virtual message_t Create(message_size_t n, some_t) override;
    virtual message_t Create(message_size_t n, mean_t) override;

   void AddTransition(message_t::size_type child, message_t::size_type parent, double weight, double u);

 protected:
    MutationModel model_;
    MutationMessageBuilder<mutation_semiring::Probability> probability_builder_;
    MutationMessageBuilder<mutation_semiring::Expectation> expectation_builder_;
    MutationMessageBuilder<mutation_semiring::Counting<0>> zero_count_builder_;
    MutationMessageBuilder<mutation_semiring::Counting<1>> one_count_builder_;
};

inline
void mutk::MutationPotential::AddTransition(message_t::size_type child, message_t::size_type parent, double weight, double u) {
    probability_builder_.AddTransition(child, parent, weight, {model_.k(), (float_t)u});
    expectation_builder_.AddTransition(child, parent, weight, {model_.k(), (float_t)u});
    zero_count_builder_.AddTransition(child, parent, weight, {model_.k(), (float_t)u});
    one_count_builder_.AddTransition(child, parent, weight, {model_.k(), (float_t)u});
}

// Possible Potentials:
// 2 x 2 -> 2 (selfing)
// 2 x 2 -> 1 (selfing)
// 1 x 1 -> 2 (selfing)
// 1 x 1 -> 1 (selfing)

// 2 x 2 -> 2
// 2 x 2 -> 1
// 2 x 1 -> 2
// 2 x 1 -> 2
// 1 x 2 -> 2
// 1 x 2 -> 1
// 1 x 1 -> 2
// 1 x 1 -> 1

// CLONING
// 2 -> 2
// 2 -> 1
// 1 -> 2
// 1 -> 1


} // namespace mutk

#endif // MUTK_POTENTIAL_HPP
