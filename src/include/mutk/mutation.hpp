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

#ifndef MUTK_MUTATION_HPP
#define MUTK_MUTATION_HPP

#include <mutk/memory.hpp>

#include <boost/container/static_vector.hpp>
#include <vector>

namespace mutk {

namespace detail {
enum struct Potential {
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

struct potential_t {
    Potential type;
    int child;
    using data_t = std::pair<int,float>;
    boost::container::static_vector<data_t, 2> parents;
    std::vector<int> shuffle;

    potential_t() = default;
    potential_t(Potential type_arg, int child_arg) : type{type_arg}, child{child_arg} {}
    potential_t(Potential type_arg, int child_arg, int par1, float dist1) :
        type{type_arg}, child{child_arg}, parents{{par1,dist1}} {}
    potential_t(Potential type_arg, int child_arg, int par1, float dist1, int par2, float dist2) :
        type{type_arg}, child{child_arg}, parents{{par1,dist1},{par2,dist2}} {}
};
} // namespace detail

namespace mutation {

// Tags
struct any_t {};
struct mean_t {};

constexpr any_t ANY = {};
constexpr mean_t MEAN = {};

class Model {
public:
    using tensor_t = mutk::Tensor<1>;
    using matrix_t = mutk::Tensor<2>;
    using potential_t = mutk::detail::potential_t;

    virtual matrix_t CreateMatrix(int n, float t, any_t) const = 0;
    virtual matrix_t CreateMatrix(int n, float t, mean_t) const = 0;
    virtual matrix_t CreateMatrix(int n, float t, int x) const = 0;

    virtual tensor_t CreatePriorDiploid(int n) const = 0;
    virtual tensor_t CreatePriorHaploid(int n) const = 0;

    mutk::Tensor<1> CreatePotential(int n, const potential_t &potential, any_t);
    mutk::Tensor<1> CreatePotential(int n, const potential_t &potential, mean_t);
    mutk::Tensor<1> CreatePotential(int n, const potential_t &potential, int x);
};

// A k-alleles model. See Lewis 2001 and Tuffley and Steel 1997.
class KAllelesModel : public Model {
public:
    using tensor_t = Model::tensor_t;
    using matrix_t = Model::matrix_t;

    KAllelesModel(double k, double theta, double hom_bias, double het_bias, double hap_bias) : 
        k_{k}, theta_{theta}, hom_bias_{hom_bias}, het_bias_{het_bias}, hap_bias_{hap_bias} {

        double e = theta/(k-1.0);

        if(k < 2.0) {
            throw std::invalid_argument("Invalid k parameter.");
        }
        if(theta < 0) {
            throw std::invalid_argument("Invalid theta parameter");
        }
        if(hom_bias > 1.0 || hom_bias < -(2.0+e)/((k-1.0)*e)) {
            throw std::invalid_argument("Invalid hom_bias parameter.");
        }
        if(het_bias > 1.0 || het_bias < -(2.0+2.0*e)/((k-2.0)*e)) {
            throw std::invalid_argument("Invalid het_bias parameter.");
        }
        if(hap_bias > 1.0 || hap_bias < -(1.0+e)/((k-1.0)*e)) {
            throw std::invalid_argument("Invalid hap_bias parameter.");
        }
    }

    virtual matrix_t CreateMatrix(int n, float t, any_t) const override;
    virtual matrix_t CreateMatrix(int n, float t, mean_t) const override;
    virtual matrix_t CreateMatrix(int n, float t, int x) const override;

    virtual tensor_t CreatePriorDiploid(int n) const override;
    virtual tensor_t CreatePriorHaploid(int n) const override;

protected:
    double k_;
    double theta_;
    double hom_bias_;
    double het_bias_;
    double hap_bias_;
};

} // namespace mutation
} // namespace mutk

#endif // MUTK_MUTATION_HPP
