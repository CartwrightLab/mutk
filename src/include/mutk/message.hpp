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

enum struct message_label_t : int {};

constexpr
message_label_t make_message_label(variable_t v, Ploidy p) {
    auto vv = static_cast<std::underlying_type_t<variable_t>>(v);
    return message_label_t{2*vv+(p == Ploidy::Diploid)};
}

constexpr
std::pair<variable_t,Ploidy> unmake_message_label(message_label_t label) {
    auto ll = static_cast<std::underlying_type_t<message_label_t>>(label);
    return std::make_pair(variable_t{ll/2}, Ploidy{(ll & 0x1)+1});
}

using float_t = float;
using message_t = xt::xarray<float_t>;
using message_shape_t = message_t::shape_type;
using message_size_t = message_t::size_type;

constexpr message_size_t num_diploids(message_size_t n) {
    return n*(n+1)/2;
}

constexpr message_size_t num_haploids(message_size_t n) {
    return n;
}

constexpr
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

constexpr
auto haploid_allele(int x) {
    return x;
}

constexpr
Ploidy message_axis_ploidy(message_label_t label) {
    return unmake_message_label(label).second;
}


constexpr message_size_t message_axis_size(message_size_t n, Ploidy h) {
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

constexpr message_size_t message_axis_size(message_size_t n, message_label_t label) {
    return message_axis_size(n, message_axis_ploidy(label));
}

// make_message(n, {pl})

} // namespace mutk

#endif // MUTK_MESSAGE_HPP
