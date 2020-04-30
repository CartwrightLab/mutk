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

#ifndef MUTK_PEELING_HPP
#define MUTK_PEELING_HPP

#include <mutk/memory.hpp>

#include <vector>
#include <utility>

#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/fill.hpp>
#include <boost/range/algorithm/fill_n.hpp>


#include <dng/genotyper.h>

#include <boost/range/adaptor/sliced.hpp>

namespace mutk {

struct workspace_t {
    GenotypeArrayVector upper; // Holds P(~Descendent_Data & G=g)
    GenotypeArrayVector lower; // Holds P( Descendent_Data | G=g)
    ParentArrayVector super; // Holds P(~Descendent_Data & G=g) for parent nodes

    bool dirty_lower = false;
    double ln_scale = 0.0;
    size_t matrix_index = 0;

    // Temporary data used by some peeling ops
    TemporaryMatrix temp_buffer;

    // Information about the pedigree
    std::size_t num_nodes = 0;
    using node_range_t = std::pair<std::size_t, std::size_t>;
    node_range_t founder_nodes,
                 germline_nodes,
                 somatic_nodes,
                 library_nodes;

    std::vector<int> ploidies;

    template<typename Rng>
    static auto make_slice(Rng& rng, node_range_t n) -> boost::sliced_range<Rng> {
        return boost::adaptors::slice(rng, n.first, n.second);
    }

    template<typename Rng>
    auto make_founder_slice(Rng& rng) -> boost::sliced_range<Rng> {
        return make_slice(rng, founder_nodes);
    }

    template<typename Rng>
    auto make_germline_slice(Rng& rng) -> boost::sliced_range<Rng> {
        return make_slice(rng, germline_nodes);
    }

    template<typename Rng>
    auto make_somatic_slice(Rng& rng) -> boost::sliced_range<Rng> {
        return make_slice(rng, somatic_nodes);
    }

    template<typename Rng>
    auto make_library_slice(Rng& rng) -> boost::sliced_range<Rng> {
        return make_slice(rng, library_nodes);
    }

    // Resize the workspace to fit a pedigree with sz nodes
    void Resize(std::size_t sz) {
        num_nodes = sz;
        upper.resize(num_nodes);
        super.resize(num_nodes);
        lower.resize(num_nodes);
        for(auto && a: lower) {
            a.setOnes(10);
        }
        temp_buffer.resize(100, 1);
        dirty_lower = false;
    }

    // Cleanup after a backwards peeling algorithm
    void Cleanup() {
        for(auto && a : lower) {
            a.setOnes();
        }
        dirty_lower = false;
    }

    // Cleanup after a backwards peeling algorithm
    // Do not update libraries since they might be set in another operation
    void CleanupFast() {
        for(std::size_t n = 0; n < somatic_nodes.second; ++n) {
            lower[n].setOnes();
        }
        dirty_lower = false;
    }

    // Set the prior probability of the founders given the reference
    void SetGermline(const GenotypeArray &prior) {
        assert(founder_nodes.first <= founder_nodes.second);
        assert(founder_nodes.second <= germline_nodes.second);

        // Set the Upper of the Founder Nodes
        for(auto i = founder_nodes.first; i < founder_nodes.second; ++i) {
            upper[i] = prior;
        }
        // Set the lowers of any germline node
        for(auto i = founder_nodes.first; i < germline_nodes.second; ++i) {
            lower[i].setOnes(prior.size());
        }
    }

    void SetGermline(const GenotypeArray &diploid_prior, const GenotypeArray &haploid_prior) {
        assert(founder_nodes.first <= founder_nodes.second);
        assert(founder_nodes.second <= germline_nodes.second);
        
        // Set the Upper of the Founder Nodes
        for(auto i = founder_nodes.first; i < founder_nodes.second; ++i) {
            assert(ploidies[i] == 2 || ploidies[i] == 1);
            if(ploidies[i] == 2) {
                upper[i] = diploid_prior;
            } else {
                upper[i] = haploid_prior;
            }
        }
        // Set the lowers of any germline node
        for(auto i = founder_nodes.first; i < germline_nodes.second; ++i) {
            assert(ploidies[i] == 2 || ploidies[i] == 1);
             if(ploidies[i] == 2) {
                lower[i].setOnes(diploid_prior.size());
            } else {
                lower[i].setOnes(haploid_prior.size());  
            }           
        }
    }

    template<typename G, typename D, typename ...A>
    void CalculateGenotypeLikelihoods(const G& gt, const D& d, A&&... args) {
        ln_scale = 0.0;
        size_t u = 0;
        for(auto pos = library_nodes.first; pos < library_nodes.second; ++pos) {
            ln_scale += gt(d[u++], std::forward<A>(args)..., ploidies[pos], &lower[pos]);
        }
    }

    void ExpGenotypeLikelihoods() {
        for(auto pos = library_nodes.first; pos < library_nodes.second; ++pos) {
            lower[pos] = lower[pos].exp();
        }
    }

    void LogGenotypeLikelihoods() {
        for(auto pos = library_nodes.first; pos < library_nodes.second; ++pos) {
            lower[pos] = lower[pos].log();
        }
    }

    // Set the genotype likelihoods into the lower values of the library_nodes.
    // Scales the genotype likelihoods as needed.
    // Input: log-likelihood values
    template<typename D>
    void SetGenotypeLikelihoods(const D& d) {
        ln_scale = 0.0;
        size_t u = 0;
        for(auto pos = library_nodes.first; pos < library_nodes.second; ++pos,++u) {
            lower[pos].resize(d[u].size());
            boost::copy(d[u], lower[pos].data());
            double temp = lower[pos].maxCoeff();
            lower[pos] = (lower[pos]-temp).exp();
            ln_scale += temp;
        }
    }    

    // Copy genotype likelihoods into the lower values of the library_nodes
    template<typename D>
    void CopyGenotypeLikelihoods(const D& d) {
        ln_scale = 0.0;
        size_t u = 0;
        for(auto pos = library_nodes.first; pos < library_nodes.second; ++pos,++u) {
            lower[pos].resize(d[u].size());
            boost::copy(d[u], lower[pos].data());
        }
    }

    // Set all genotype likelihoods to 1
    void ClearGenotypeLikelihoods(size_t num_of_obs_alleles) {
        ln_scale = 0.0;
        size_t hap_sz = num_of_obs_alleles;
        size_t gt_sz = hap_sz*(hap_sz+1)/2;
        for(auto pos = library_nodes.first; pos < library_nodes.second; ++pos) {
            assert(ploidies[pos] == 1 || ploidies[pos] == 2);
            lower[pos].setOnes( (ploidies[pos] == 2) ? gt_sz : hap_sz );
        }
    }   
};

/* Peeling Operators
PeelUp:       low[o] = mat %*% low[a]
PeelUpStar:   low[o] = (mat %*% low[a])*low[b]
PeelDown:     upp[o] = mat  %*% upp[a]
PeelStar:     upp[o] = upp[a] * low[b]
PeelMate:     low[o] = upp[a] %*% reshaped[low[b]]
*/


class PeelOp {
public:
    PeelOp(int o, int a, int b, TransitionMatrix mat) :
        output_{o}, a_{a}, b_{b}, mat_{std::move(mat)}
    {}

    virtual void operator()(workspace_t *work) const = 0;

    // Node indexes
    int output_;
    int a_;
    int b_;

    TransitionMatrix mat_;
};

class PeelUp : public PeelOp {
public:
    PeelUp(int o, int a, TransitionMatrix mat) :
    PeelOp{o, a, -1, std::move(mat)} {}
    virtual void operator()(workspace_t *work) const final;
};

class PeelUpStar : public PeelOp {
public:
    using PeelOp::PeelOp;
    virtual void operator()(workspace_t *work) const final;
};

class PeelDown : public PeelOp {
public:
    PeelDown(int o, int a, int b, TransitionMatrix mat) :
    PeelOp{o, a, b, std::move(mat)} {
        mat.transposeInPlace();
    }
    virtual void operator()(workspace_t *work) const final;
};

void to_father(workspace_t &work, const family_members_t &family,
               const TransitionMatrixVector &mat);
void to_mother(workspace_t &work, const family_members_t &family,
               const TransitionMatrixVector &mat);
void to_child(workspace_t &work, const family_members_t &family,
              const TransitionMatrixVector &mat);

// Fast versions that can be used on "dirty" workspaces
void up_fast(workspace_t &work, const family_members_t &family,
             const TransitionMatrixVector &mat);
void down_fast(workspace_t &work, const family_members_t &family,
               const TransitionMatrixVector &mat);
void to_father_fast(workspace_t &work, const family_members_t &family,
                    const TransitionMatrixVector &mat);
void to_mother_fast(workspace_t &work, const family_members_t &family,
                    const TransitionMatrixVector &mat);
void to_child_fast(workspace_t &work, const family_members_t &family,
                   const TransitionMatrixVector &mat);

// Reverse versions that can be used for backwards algorithms
void up_reverse(workspace_t &work, const family_members_t &family,
                const TransitionMatrixVector &mat);
void down_reverse(workspace_t &work, const family_members_t &family,
                  const TransitionMatrixVector &mat);
void to_father_reverse(workspace_t &work, const family_members_t &family,
                       const TransitionMatrixVector &mat);
void to_mother_reverse(workspace_t &work, const family_members_t &family,
                       const TransitionMatrixVector &mat);
void to_child_reverse(workspace_t &work, const family_members_t &family,
                      const TransitionMatrixVector &mat);

typedef decltype(&down) function_t;

struct info_t {
    bool writes_lower;
    int writes_to;
};

enum struct Op { 
    UP=0, DOWN, TOFATHER, TOMOTHER, TOCHILD,
    UPFAST, DOWNFAST, TOFATHERFAST, TOMOTHERFAST,
    TOCHILDFAST,
    NUM // Total number of possible forward operations
};

constexpr info_t info[(int)Op::NUM] = {
    /* Up           */ {true,  0},
    /* Down         */ {false, 1},
    /* ToFather     */ {true,  0},
    /* ToMother     */ {true,  1},
    /* ToChild      */ {false, 2},
    /* UpFast       */ {true,  0},
    /* DownFast     */ {false, 1},
    /* ToFatherFast */ {true,  0},
    /* ToMotherFast */ {true,  1},
    /* ToChildFast  */ {false, 2}
};

// TODO: Write test case to check that peeling ops are in the right order.
constexpr function_t functions[(int)Op::NUM] = {
    &up, &down, &to_father, &to_mother, &to_child,
    &up_fast, &down_fast, &to_father_fast, &to_mother_fast,
    &to_child_fast
};

constexpr function_t reverse_functions[(int)Op::NUM] = {
    &up_reverse, &down_reverse, &to_father_reverse,
    &to_mother_reverse, &to_child_reverse,
    &up_reverse, &down_reverse, &to_father_reverse,
    &to_mother_reverse, &to_child_reverse
};

} // namespace dng::peel
} // namespace dng

#endif // DNG_PEELING_H
