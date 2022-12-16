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

#include <vector>

#include "message.hpp"
#include "mutation.hpp"

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
    using labels_t = message_labels_t;

    Potential(labels_t labels) : labels_{std::move(labels)} {}

    virtual ~Potential() = default;

    enum struct any_t  : int {};
    enum struct mean_t : int {};
    enum struct some_t : int {};

    virtual message_t Create(message_size_t n, any_t) = 0;
    virtual message_t Create(message_size_t n, some_t) = 0;
    virtual message_t Create(message_size_t n, mean_t) = 0;

    static constexpr any_t ANY{0};
    static constexpr mean_t MEAN{1};
    static constexpr some_t ZERO{0};
    static constexpr some_t ONE{1};
    static constexpr some_t TWO{2};

 protected:
    labels_t labels_;
};

inline
constexpr auto operator+(Potential::some_t value) {
    return static_cast<std::underlying_type_t<Potential::some_t>>(value);
}

class UnitPotential : public Potential {
 public:
    UnitPotential(labels_t labels) : Potential(std::move(labels)) { }

    virtual message_t Create(message_size_t n, any_t) override;
    virtual message_t Create(message_size_t n, some_t) override;
    virtual message_t Create(message_size_t n, mean_t) override;
};

class MutationPotential : public Potential {
 public:
    MutationPotential(labels_t labels, const MutationModel & model) : Potential(std::move(labels)),
        model_(model) { }

 protected:
    MutationModel model_;
};

template<Ploidy P, Ploidy C>
class CloningPotential : public MutationPotential {
 public:
    CloningPotential(labels_t labels, const MutationModel & model,
        float_t u) : MutationPotential(std::move(labels), model),
        u_(u) { }

    virtual message_t Create(message_size_t n, any_t val) override;
    virtual message_t Create(message_size_t n, some_t) override;
    virtual message_t Create(message_size_t n, mean_t) override;    

 protected:
    template<class Arg>
    message_t DoCreate(message_size_t n, Arg);

    inline auto MessageShape(message_size_t n) {
        auto d1 = message_dimension_length(n, P);
        auto d2 = message_dimension_length(n, C);
        return std::array<message_size_t, 2>({d1, d2});
    }

    float_t u_;
};

class SelfingPotential : public MutationPotential {
 public:
    SelfingPotential(labels_t labels, const MutationModel & model,
        float_t len_a, float_t len_b) : MutationPotential(std::move(labels), model),
        len_a_(len_a), len_b_(len_a) { }

    virtual message_t Create(message_size_t n, any_t) override;
    virtual message_t Create(message_size_t n, some_t) override;
    virtual message_t Create(message_size_t n, mean_t) override;    

 protected:
    float_t len_a_;
    float_t len_b_;
};

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

template<Ploidy A, Ploidy B>
message_t CloningPotential<A,B>::Create(size_t n, any_t a) {
    return DoCreate(n, a);
}

template<Ploidy A, Ploidy B>
message_t CloningPotential<A,B>::Create(size_t n, mean_t a) {
    return DoCreate(n, a);
}

template<Ploidy A, Ploidy B>
message_t CloningPotential<A,B>::Create(size_t n, some_t a) {
    return DoCreate(n, a);
}

} // namespace mutk

#endif // MUTK_POTENTIAL_HPP
