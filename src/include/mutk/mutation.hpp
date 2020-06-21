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
    using tensor_t = mutk::Tensor<2>;
    using potential_t = mutk::detail::potential_t;

    virtual tensor_t CreateMatrix(int n, float t, any_t) const = 0;
    virtual tensor_t CreateMatrix(int n, float t, mean_t) const = 0;
    virtual tensor_t CreateMatrix(int n, float t, int x) const = 0;

    mutk::Tensor<1> CreatePotential(int n, const potential_t &potential, any_t);
};

// A k-alleles model. See Lewis 2001 and Tuffley and Steel 1997.
class KAllelesModel : public Model {
public:
    using tensor_t = Model::tensor_t;

    KAllelesModel(double k) : k_{k} {
        assert(k_ >= 2.0);        
    }

    virtual tensor_t CreateMatrix(int n, float t, any_t) const override;
    virtual tensor_t CreateMatrix(int n, float t, mean_t) const override;
    virtual tensor_t CreateMatrix(int n, float t, int x) const override;

protected:
    double k_;
};



/*
// k-alleles model from Watterson and Guess (1977) https://doi.org/10.1016/0040-5809(77)90023-5

inline
tensor_t population_prior_diploid(int num_obs_alleles, double theta, double hom_bias, double het_bias,
    double kalleles) {
    assert(num_obs_alleles >= 0);

    double k = kalleles;
    double e = theta/(k-1.0);

    double p_hom = (1.0+e)/(1.0+k*e);
    double p_hetk = e/(1.0+k*e);

    double p_RR = p_hom*(2.0+e+(k-1.0)*e*hom_bias)/(2.0+k*e);
    double p_AA = p_hom*(e-e*hom_bias)/(2.0+k*e);

    double p_RA = p_hetk*(2.0+2.0*e+(k-2.0)*e*het_bias)/(2.0+k*e);
    double p_AB = p_hetk*(2.0*e-2.0*e*het_bias)/(2.0+k*e);

    tensor_t ret{num_obs_alleles*(num_obs_alleles+1)/2};

    int n=0;
    for(int i=0;i<num_obs_alleles;++i) {
        for(int j=0;j<i;++j) {
            ret(n++) = (j==0 || i==0) ? p_RA : p_AB;
        }
        ret(n++) = (i==0) ? p_RR : p_AA;
    }

    return ret;
}

inline
tensor_t population_prior_haploid(int num_obs_alleles, double theta, double hap_bias,
    double kalleles) {
    assert(num_obs_alleles >= 1);

    double k = kalleles;
    double e = theta/(k-1.0);

    double p_R = (1.0+e+(k-1.0)*e*hap_bias)/(1.0+k*e);
    double p_A = (e-e*hap_bias)/(1.0+k*e);

    tensor_t ret{num_obs_alleles};
    ret(0) = p_R;
    for(int n=1;n<num_obs_alleles;++n) {
        ret(n) = p_A;
    }

    return ret;
}

inline bool population_prior_check(double theta, double hom_bias, double het_bias, double hap_bias, double kalleles) {
    double k = kalleles;
    double e = theta/(k-1.0);

    if(e < 0) {
        return false;
    }
    if(k < 2.0) {
        return false;
    }
    if(hom_bias > 1.0 || hom_bias < -(2.0+e)/((k-1.0)*e)) {
        return false;
    }
    if(het_bias > 1.0 || het_bias < -(2.0+2.0*e)/((k-2.0)*e)) {
        return false;
    }
    if(hap_bias > 1.0 || hap_bias < -(1.0+e)/((k-1.0)*e)) {
        return false;
    }
    return true;
}

*/

} // namespace mutation
} // namespace mutk

#endif // MUTK_MUTATION_HPP
