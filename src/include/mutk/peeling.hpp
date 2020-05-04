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

#include <boost/range/adaptor/sliced.hpp>

namespace mutk {

struct workspace_t;

/* Peeling Operators
PeelUp:       low[o] = mat %*% low[a]
PeelUpStar:   low[o] = (mat %*% low[a])*low[b]
PeelDown:     upp[o] = mat  %*% upp[a]
PeelStar:     upp[o] = upp[a] * low[b]
PeelMate:     low[o] = reshaped[low[a]] %*% upp[b]
*/

class PeelOp {
public:
    PeelOp(int o, int a, int b) :
        output_{o}, a_{a}, b_{b}
    {}

    virtual void operator()(workspace_t *work, const TransitionMatrix &mat) const = 0;

    // Node indexes
    int output_;
    int a_;
    int b_;
};

class PeelUp : public PeelOp {
public:
    PeelUp(int o, int a) :
    PeelOp{o, a, -1} {}
    virtual void operator()(workspace_t *work, const TransitionMatrix &mat) const final;
};

class PeelUpStar : public PeelOp {
public:
    using PeelOp::PeelOp;
    virtual void operator()(workspace_t *work, const TransitionMatrix &mat) const final;
};

class PeelDown : public PeelOp {
public:
    using PeelOp::PeelOp;
    virtual void operator()(workspace_t *work, const TransitionMatrix &mat) const final;
};

class PeelStar : public PeelOp {
public:
    using PeelOp::PeelOp;
    virtual void operator()(workspace_t *work, const TransitionMatrix &mat) const final;
};

class PeelMate : public PeelOp {
public:
    using PeelOp::PeelOp;
    virtual void operator()(workspace_t *work, const TransitionMatrix &mat) const final;
};


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

} // namespace MUTK

#endif // MUTK_PEELING_HPP
