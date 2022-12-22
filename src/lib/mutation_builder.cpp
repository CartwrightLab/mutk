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

#include "unit_testing.hpp"
#include "mutation_testing.hpp"

#include <cmath>
#include <iterator>
#include <utility>
#include <mutk/message.hpp>

#include <boost/math/special_functions/binomial.hpp>

class MutationSemiring {
public:
    MutationSemiring(double k, double u) : k_(k), u_(u) {
        double beta = k_/(k_-1.0);
        pij_ = -1.0/k_*expm1(-beta*u_);
        pii_ = exp(-beta*u_) + pij_;
    }

    using value_type = double;

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

    value_type operator()(int a, int b, double w) const {
        return w * ((a==b) ? pii_ : pij_);
    }

private:
    double k_;
    double u_;
    double pii_;
    double pij_;
};

// Temporary Name
template<class T>
class MutationBuilder { 
public:
    using mutation_type = T;
    using message_type = mutk::message_t;

    explicit MutationBuilder(std::vector<int> ploidies) : ploidies_{std::move(ploidies)} {
        if(!ploidies_.empty()) {
            child_ploidy_ = ploidies_[0];
            parents_ploidy_ = std::accumulate(std::next(ploidies_.begin()),
                ploidies_.end(), parents_ploidy_);
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


// Find right most minimum number. Increase that number by one. Set all values to the left of it to zero

template<class BidirIt, class ValueType>
bool next_multiset(BidirIt first, BidirIt last, ValueType n) {
    if(first == last) {
        return false;
    }
    auto it = std::find_if(first, last, [&](auto cit){ return *first < cit; });
    std::advance(it, -1);
    if(++(*it) == n) {
        std::fill(first, last, 0);
        return false;
    }
    std::fill(first, it, 0);
    return true;
}

template<class T>
auto MutationBuilder<T>::Create(int n) const -> message_type {
    message_type::shape_type shape = Shape(n);
    message_type msg(shape);

    using BidirIt = typename std::vector<int>::iterator;

    std::vector<int> coords(child_ploidy_+parents_ploidy_, 0);

    std::vector<std::pair<BidirIt,BidirIt>> partitions;
    int pos = 0;
    for(int p : ploidies_) {
        partitions.emplace_back(coords.begin()+pos, coords.begin()+pos+p);
        pos += p;
    }

    auto do_permute = [&](auto &val) {
        return std::next_permutation(val.first, val.second);
    };
    auto do_next_multiset = [&](auto &val) {
        return next_multiset(val.first, val.second, n);
    };

    message_type::size_type idx = 0;
    do {
        auto total = mutation_type::Zero();
        int counter = 0;
        do {
            do {
                auto temp = mutation_type::One();
                for(int x = 0; x < child_ploidy_; ++x) {
                    auto value = mutation_type::Zero();
                    for(auto &&par : transitions_[x]) {
                        value = mutation_type::Plus(value, par.mu(coords[x], coords[par.parent], par.weight));
                    }
                    temp = mutation_type::Times(temp, value);
                }
                total = mutation_type::Plus(total, temp);
            } while(do_permute(partitions[0]));
            counter += 1;
        } while(std::any_of(std::next(partitions.begin()), partitions.end(), do_permute));

        msg.flat(idx++) = total / counter;
    } while(std::any_of(partitions.begin(), partitions.end(), do_next_multiset));

    return msg;
}

template<class T>
auto MutationBuilder<T>::Shape(int n) const -> message_type::shape_type {
    message_type::shape_type ret;
    for(auto v : ploidies_) {
        ret.push_back(boost::math::binomial_coefficient<double>(n+v-1, v));
    }
    std::reverse(ret.begin(), ret.end());
    return ret;
}

TEST_CASE("MutationBuilder") {
    using Builder = MutationBuilder<MutationSemiring>;

    using S = mutk::message_t::shape_type;

    {
        Builder builder({2,2});

        builder.AddTransition(0, 2, 1.0, MutationSemiring(4,0.001));
        builder.AddTransition(1, 3, 1.0, MutationSemiring(4,0.001));

        double k = 4;
        double u = 0.001;
        double beta = k/(k-1.0);
        double pij = -1.0/k*expm1(-beta*u);
        double pii = exp(-beta*u) + pij;

        auto expected = mutk::message_t::from_shape({10,10});

        for(int a=0,i=0;a<4;++a) {
            for(int b=0;b<=a;++b,++i) {
                for(int x=0,j=0;x<4;++x) {
                    for(int y=0;y<=x;++y,++j) {
                        // ab -> xy
                        double ax = (a == x) ? pii : pij;
                        double by = (b == y) ? pii : pij;
                        double ay = (a == y) ? pii : pij;
                        double bx = (b == x) ? pii : pij;
                        double ret = ax*by;
                        if(x != y) {
                            ret += ay*bx;
                        }
                        expected(i,j) = ret;
                    }
                }
            }
        }

        auto msg = builder.Create(4);
        CHECK_EQ_RANGES(msg, expected);
    }
    {
        Builder builder({0,2,2});

        auto expected = mutk::message_t::from_shape({10,10});
        expected.fill(1.0);
        auto msg = builder.Create(4);
        CHECK_EQ_RANGES(msg, expected);
    }
    {
        Builder builder({2,2,1});

        builder.AddTransition(0, 2, 0.5, MutationSemiring(4,0.0001));
        builder.AddTransition(0, 3, 0.5, MutationSemiring(4,0.0001));
        builder.AddTransition(1, 4, 1.0, MutationSemiring(4,0.001));

        auto msg = builder.Create(2);

        double k = 4;
        double u = 0.0001;
        double v = 0.001;
        double beta = k/(k-1.0);
        double piju = -1.0/k*expm1(-beta*u);
        double piiu = exp(-beta*u) + piju;
        double pijv = -1.0/k*expm1(-beta*v);
        double piiv = exp(-beta*v) + pijv;

        auto expected = mutk::message_t::from_shape({2,3,3});

        int n = 2;
        for(int a=0,i=0;a<n;++a) {
            for(int b=0;b<=a;++b,++i) {
                for(int c=0,j=0;c<n;++c,++j) {
                    for(int x=0,k=0;x<n;++x) {
                        for(int y=0;y<=x;++y,++k) {
                            // ab -> xy
                            double ax = (a == x) ? piiu : piju;
                            double bx = (b == x) ? piiu : piju;
                            double abx = 0.5*(ax+bx);
                            double cy = (c == y) ? piiv : pijv;
                            double ret = abx*cy;
                            if(x != y) {
                                ax = (a == y) ? piiu : piju;
                                bx = (b == y) ? piiu : piju;
                                abx = 0.5*(ax+bx);
                                cy = (c == x) ? piiv : pijv;
                                ret += abx*cy;
                            }
                            expected(j,i,k) = ret;
                        }
                    }
                }
            }
        }

        CHECK_EQ_RANGES(msg, expected);
    }

    //std::cout << msg << std::endl;
}
