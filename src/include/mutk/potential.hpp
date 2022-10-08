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
#include <boost/container/static_vector.hpp>

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

// A strongly-type individual id.
// Use unitary + to do a static cast.
enum struct IndyId : int {};
constexpr auto operator+(IndyId id) {
    return static_cast<std::underlying_type_t<IndyId>>(id);
}

struct potential_t {
    PotentialType type;
    IndyId child;
    using data_t = std::pair<IndyId,float>;
    boost::container::static_vector<data_t, 2> parents;
    std::vector<int> axes;

    potential_t() = default;
    potential_t(PotentialType type_arg, IndyId child_arg) : type{type_arg}, child{child_arg} {}
    potential_t(PotentialType type_arg, IndyId child_arg, IndyId par1, float dist1) :
        type{type_arg}, child{child_arg}, parents{{par1,dist1}} {}
    potential_t(PotentialType type_arg, IndyId child_arg, IndyId par1, float dist1, IndyId par2, float dist2) :
        type{type_arg}, child{child_arg}, parents{{par1,dist1},{par2,dist2}} {}
};

} // namespace mutk

#endif // MUTK_POTENTIAL_HPP
