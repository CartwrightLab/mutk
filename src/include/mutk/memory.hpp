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

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xgenerator.hpp>

namespace mutk {

using float_t = float;
using tensor_t = xt::xarray<float_t>;
using shape_t = tensor_t::shape_type;
using strides_t = tensor_t::strides_type;

inline
constexpr tensor_t::size_type num_diploids(tensor_t::size_type n) {
    return n*(n+1)/2;
}

template<int N>
constexpr tensor_t::size_type dim_width(tensor_t::size_type n) {
    static_assert(N == 1 || N == 2);
    if constexpr(N == 1) {
        return n;
    } else {
        return num_diploids(n);
    }
}

// namespace detail {
//     // helper function copied from detail namespace of xtensor library
//     template <class Functor, class S>
//     inline auto make_xgenerator(Functor&& f, S&& shape) noexcept {
//         using type = xt::xgenerator<Functor, std::invoke_result_t<Functor, S>, std::decay_t<S>>;
//         return type(std::forward<Functor>(f), std::forward<S>(shape));
//     }
// }

// template<std::size_t N>
// using Tensor = Eigen::Tensor<float,N>;

// template<typename Arg1, typename ...Args>
// typename Eigen::array<Eigen::DenseIndex,1+sizeof...(Args)>
// tensor_dims(Arg1 arg1, Args... args) {
//     return {arg1,args...};
// }

// template<typename T, typename ...Args>
// auto wrap_tensor(T *p, Args... args) {
//     return Eigen::TensorMap<Eigen::Tensor<T,sizeof...(Args)>>(p, args...);
// }

// template<typename T, std::size_t N>
// auto wrap_tensor(T *p, Eigen::array<Eigen::DenseIndex,N> arg) {
//     return Eigen::TensorMap<Eigen::Tensor<T,N>>(p, arg);
// }

// template<typename A, typename B>
// inline auto kronecker_product_coef(const A &a, const B &b, std::size_t i,
//                               std::size_t j) {
//     assert(i < a.rows()*b.rows() && j < a.cols()*b.cols());
//     return a(i / b.rows(), j / b.cols()) * b(i % b.rows(), j % b.cols());
// }

// From http://stackoverflow.com/a/16287999
// template <typename T, int R, int C>
// inline T sum_kahan(const Eigen::Array<T, R, C> &xs) {
//     if(xs.size() == 0) {
//         return 0;
//     }
//     T sumP(0), sumN(0), tP(0), tN(0);
//     T cP(0), cN(0), yP(0), yN(0);
//     for(size_t i = 0, size = xs.size(); i < size; i++) {
//         T temporary = (*(xs.data() + i));
//         if(temporary > 0) {
//             yP = temporary - cP;
//             tP = sumP + yP;
//             cP = (tP - sumP) - yP;
//             sumP = tP;
//         } else {
//             yN = temporary - cN;
//             tN = sumN + yN;
//             cN = (tN - sumN) - yN;
//             sumN = tN;
//         }
//     }
//     return sumP + sumN;
// }

} // namespace mutk

// for debugging purposes
namespace std {
inline
std::ostream& operator<< (std::ostream& os, const mutk::shape_t& value) {
    os << "{";
    for(size_t i = 0; i < value.size(); ++i) {
        if(i > 0) {
            os << ", ";
        }
        os << value[i];
    }
    os << "}";
    return os;
}
}



#endif // MUTK_MEMORY_HPP
