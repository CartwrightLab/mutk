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

#ifndef MUTK_MEMORY_HPP
#define MUTK_MEMORY_HPP

#include <cstdint>
#include <cfloat>
#include <algorithm>
#include <initializer_list>

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/CXX11/Tensor>

namespace mutk {

template<std::size_t N>
using Tensor = Eigen::Tensor<float,N>;




template<typename A, typename B>
inline auto kronecker_product_coef(const A &a, const B &b, std::size_t i,
                              std::size_t j) {
    assert(i < a.rows()*b.rows() && j < a.cols()*b.cols());
    return a(i / b.rows(), j / b.cols()) * b(i % b.rows(), j % b.cols());
}

// From http://stackoverflow.com/a/16287999
template <typename T, int R, int C>
inline T sum_kahan(const Eigen::Array<T, R, C> &xs) {
    if(xs.size() == 0) {
        return 0;
    }
    T sumP(0), sumN(0), tP(0), tN(0);
    T cP(0), cN(0), yP(0), yN(0);
    for(size_t i = 0, size = xs.size(); i < size; i++) {
        T temporary = (*(xs.data() + i));
        if(temporary > 0) {
            yP = temporary - cP;
            tP = sumP + yP;
            cP = (tP - sumP) - yP;
            sumP = tP;
        } else {
            yN = temporary - cN;
            tN = sumN + yN;
            cN = (tN - sumN) - yN;
            sumN = tN;
        }
    }
    return sumP + sumN;
}

} // namespace mutk

#endif // MUTK_MEMORY_HPP
