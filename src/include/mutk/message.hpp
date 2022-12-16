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

#ifndef MUTK_MESSAGE_HPP
#define MUTK_MESSAGE_HPP

#include <utility>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xgenerator.hpp>

#include <boost/container/flat_map.hpp>

namespace mutk {

enum struct variable_t : int {};
constexpr auto operator+(variable_t value) {
    return static_cast<std::underlying_type_t<variable_t>>(value);
}

enum struct Ploidy : int {
    Noploid = 0,
    Haploid = 1,
    Diploid = 2
};

using float_t = float;
using message_t = xt::xarray<float_t>;
using message_shape_t = message_t::shape_type;
using message_size_t = message_t::size_type;

using message_labels_t = boost::container::flat_map<variable_t, Ploidy>;

inline
constexpr message_size_t num_diploids(message_size_t n) {
    return n*(n+1)/2;
}

inline
constexpr message_size_t num_haploids(message_size_t n) {
    return n;
}

inline constexpr
auto diploid_alleles(int x) {
    constexpr int ALLELE[][2] = {
        {0,0},
        {0,1}, {1,1},
        {0,2}, {1,2}, {2,2},
        {0,3}, {1,3}, {2,3}, {3,3},
        {0,4}, {1,4}, {2,4}, {3,4}, {4,4},
        {0,5}, {1,5}, {2,5}, {3,5}, {4,5}, {5,5},
    };
    return std::make_pair(ALLELE[x][0], ALLELE[x][1]);
}

inline constexpr
auto haploid_allele(int x) {
    return x;
}

inline
constexpr message_size_t message_dimension_length(message_size_t n, Ploidy h) {
    switch(h) {
     case Ploidy::Noploid:
        return 0;
     case Ploidy::Haploid:
        return num_haploids(n);
     case Ploidy::Diploid:
        return num_diploids(n);
    };
    return 0;
}

inline
message_shape_t message_shape(message_size_t n, const message_labels_t & labels) {
    const auto & seq = labels.sequence();
    message_shape_t ret;
    ret.resize(seq.size());
    for(message_size_t i=0; i < seq.size(); ++i) {
        ret[i] = message_dimension_length(n, seq[i].second);
    }
    return ret;
}

// make_message(n, {pl})

} // namespace mutk

#endif // MUTK_MESSAGE_HPP
