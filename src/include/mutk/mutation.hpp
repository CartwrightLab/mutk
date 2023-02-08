/*
# Copyright (c) 2022 Reed A. Cartwright <racartwright@gmail.com>
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

#include <stdexcept>

#include <cmath>
#include <iterator>
#include <utility>

#include <boost/math/special_functions/binomial.hpp>

#include "message.hpp"

namespace mutk {

// A k-alleles model. See Lewis 2001 and Tuffley and Steel 1997.
//
// P(j|i) = 1/k * (1-exp(-beta*t)) + I(j == i)*exp(-beta*t)
//    where beta = k/(k-1)
//    expected number of mutations = t
//    expected number of events = beta*t
//
//    1/k mutations are virtual and (k-1)/k are real
//
//    Given at least 1 event, the expected number of events is 
//      (beta*t)/(1-exp(-beta*t))
//
//    And expected number of mutations is t/(1-exp(-beta*t))
//
class MutationModel {
public:
    using array_t = mutk::message_t;

    MutationModel(float_t k, float_t theta, float_t hom_bias, float_t het_bias, float_t hap_bias) : 
        k_{k}, theta_{theta}, hom_bias_{hom_bias}, het_bias_{het_bias}, hap_bias_{hap_bias} {

        float_t e = theta/(k-1.0f);

        if(k < 2.0f) {
            throw std::invalid_argument("Invalid k parameter.");
        }
        if(theta < 0.0f) {
            throw std::invalid_argument("Invalid theta parameter");
        }
        if(hom_bias > 1.0f || hom_bias < -(2.0f+e)/((k-1.0f)*e)) {
            throw std::invalid_argument("Invalid hom_bias parameter.");
        }
        if(het_bias > 1.0f || het_bias < -(2.0f+2.0f*e)/((k-2.0f)*e)) {
            throw std::invalid_argument("Invalid het_bias parameter.");
        }
        if(hap_bias > 1.0f || hap_bias < -(1.0f+e)/((k-1.0f)*e)) {
            throw std::invalid_argument("Invalid hap_bias parameter.");
        }
    }

    float_t k() const { return k_; }
    float_t theta() const { return theta_; }
    float_t hom_bias() const { return hom_bias_; }
    float_t het_bias() const { return het_bias_; }
    float_t hap_bias() const { return hap_bias_; }

protected:
    float_t k_;
    float_t theta_;
    float_t hom_bias_;
    float_t het_bias_;
    float_t hap_bias_;
};

namespace mutation_semiring {

class Probability {
public:
    Probability(float_t k, float_t u) : k_(k), u_(u) {
        double beta = k_/(k_-1.0);
        double p = -expm1(-beta*u_);
        pij_ = 1.0/k_*p;
        pii_ = exp(-beta*u_) + pij_;
    }

    using value_type = float_t;

    value_type operator()(int a, int b, float_t w) const {
        return w * ((a==b) ? pii_ : pij_);
    }

    static value_type One() { return 1.0; }
    static value_type Zero() { return 0.0; }

    // Combination
    static value_type Plus(value_type lhs, value_type rhs) {
        return lhs + rhs;
    }

    // Aggregation
    static value_type Times(value_type lhs, value_type rhs) {
        return lhs * rhs;
    }

private:
    float_t k_;
    float_t u_;
    float_t pii_;
    float_t pij_;
};

// https://www.cs.jhu.edu/~jason/papers/li+eisner.emnlp09.pdf
// Value type is pair <p,r>
//     p = P(j|i)
//     r = E[num of mutations | i,j]*P(j|i)
class Expectation {
public:
    Expectation(float_t k, float_t u) : k_(k), u_(u) {
        double beta = k_/(k_-1.0);
        double p = -expm1(-beta*u_);

        pii_ = exp(-beta*u_) + 1.0/k_*p;
        rii_ = p * u_/k_;

        pij_ = 1.0/k_*p;
        rij_ = (k_-p)/(k_-1.0) * u_/k_;
    }

    using value_type = std::array<float_t,2>;

    value_type operator()(int a, int b, float_t w) const {
        return ((a==b) ? value_type{w * pii_, w * rii_} :
                         value_type{w * pij_, w * rij_});
    }

    static value_type One() { return {1.0,0.0}; }
    static value_type Zero() { return {0.0,0.0}; }

    // Combination
    static value_type Plus(value_type lhs, value_type rhs) {
        return {lhs[0] + rhs[0], lhs[1] + rhs[1]};
    }

    // Aggregation
    static value_type Times(value_type lhs, value_type rhs) {
        return {lhs[0]*rhs[0], lhs[0]*rhs[1] + rhs[0]*lhs[1]};
    }

private:
    float_t k_;
    float_t u_;
    float_t pii_;
    float_t rii_;
    float_t pij_;
    float_t rij_;
};

} // namespace mutation_semiring

template<class T>
class MutationMessageBuilder { 
public:
    using mutation_type = T;
    using message_type = mutk::message_t;

    explicit MutationMessageBuilder(std::vector<int> ploidies) : ploidies_{std::move(ploidies)} {
        if(!ploidies_.empty()) {
            child_ploidy_ = ploidies_.back();
            parents_ploidy_ = std::accumulate(ploidies_.begin(),
                ploidies_.end(), parents_ploidy_) - child_ploidy_;
            transitions_.resize(child_ploidy_);
        }
    }

    void AddTransition(message_type::size_type child, message_type::size_type parent, double weight, mutation_type mu) {
        transitions_[child].push_back({weight, parent, mu});
    }

    message_type::shape_type Shape(int n) const;

    message_type Create(int n) const;

private:
    int child_ploidy_{0};
    int parents_ploidy_{0};
    std::vector<int> ploidies_;

    struct transition_t {
        double weight;
        message_type::size_type parent;
        mutation_type mu;
    };

    std::vector<std::vector<transition_t>> transitions_;
};

// If k goes from [0 to N)
// The index for genotype k1 <= k2 <= ... <= kP is
// Sum_{m=1}^P (kM+m-1) choose m

// The number of possible genotypes with ploidy P and alleles N is
// (N+P-1) choose P

namespace detail {
// Find right most minimum number. Increase that number by one. Set all values to the left of it to zero
template<class BidirIt, class ValueType>
bool next_multiset(BidirIt first, BidirIt last) {
    if(first == last) {
        return false;
    }
    auto it = std::find_if(first, last, [&](auto cit){ return *first < cit; });
    std::advance(it, -1);
    ++(*it);
    std::fill(first, it, 0);
    return true;
}

}

template<class T>
auto MutationMessageBuilder<T>::Create(int n) const -> message_type {
    message_type::shape_type shape = Shape(n);
    message_type msg(shape);

    using BidirIt = typename std::vector<int>::iterator;

    std::vector<int> coords(child_ploidy_+parents_ploidy_, 0);

    std::vector<std::pair<BidirIt,BidirIt>> partitions;
    int pos = parents_ploidy_+child_ploidy_;
    for(auto it = ploidies_.rbegin(); it != ploidies_.rend(); ++it) {
        pos -= *it;
        partitions.emplace_back(coords.begin()+pos, coords.begin()+pos+*it);
    }

    auto do_next_order = [&](auto &val) {
        return std::next_permutation(val.first, val.second);
    };
    auto do_next_multiset = [&](auto &val) {
        // The multiset iteration needs to reset to all 0s once the last element is equal to n.
        // Return false when it resets to signal a new round is starting.
        if(detail::next_multiset(val.first, val.second)) {
            if(*std::prev(val.second) < n) {
                return true;
            }
            std::fill(val.first, val.second, 0);
        }
        return false;
    };

    message_type::size_type idx = 0;
    do {
        auto total = mutation_type::Zero();
        int counter = 0;
        do {
            // Loop through phased parent genotypes
            do {
                // Loop through phased child genotypes
                auto temp = mutation_type::One();
                for(int x = 0; x < child_ploidy_; ++x) {
                    auto value = mutation_type::Zero();
                    for(auto &&par : transitions_[x]) {
                        value = mutation_type::Plus(value, par.mu(coords[x+parents_ploidy_],
                            coords[par.parent], par.weight));
                    }
                    temp = mutation_type::Times(temp, value);
                }
                total = mutation_type::Plus(total, temp);
            } while(do_next_order(partitions[0]));
            counter += 1;
        } while(std::any_of(std::next(partitions.begin()), partitions.end(), do_next_order));

        msg.flat(idx++) = total / counter;
    } while(std::any_of(partitions.begin(), partitions.end(), do_next_multiset));

    return msg;
}

template<class T>
auto MutationMessageBuilder<T>::Shape(int n) const -> message_type::shape_type {
    message_type::shape_type ret;
    for(auto v : ploidies_) {
        if(v > 0) {
            ret.push_back(boost::math::binomial_coefficient<double>(n+v-1, v));
        }
    }
    return ret;
}

} // namespace mutk

#endif // MUTK_MUTATION_HPP
