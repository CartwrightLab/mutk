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

#include <mutk/potential.hpp>

using mutk::message_t;
using mutk::Ploidy;

// Use some template magic to allow maximum deduplication of code.
struct mutk::MatingPotential::Impl {
    Impl(const mutk::MatingPotential &p) : pot{p} {}

    template<class Arg>
    inline
    message_t operator()(size_t n, Arg a);

    // Use nested class to mimic partial template specialization for functions
    template<int N>
    struct Create {
        template<class Arg>
        static message_t call(const Impl &impl, size_t n, Arg a);
    };

    inline
    auto CreateModelMatrix(size_t n, mutk::Potential::any_t, float t) const {
        return pot.model_.CreateTransitionMatrix(n, t);
    }

    inline
    auto CreateModelMatrix(size_t n, mutk::Potential::mean_t a, float t) const {
        auto aa = static_cast<std::underlying_type_t<mutk::Potential::mean_t>>(a);
        return (aa > 0) ? pot.model_.CreateMeanMatrix(n, t) :
                          pot.model_.CreateTransitionMatrix(n, t);
    }

    inline
    auto CreateModelMatrix(size_t n, mutk::Potential::some_t a, float t) const {
        auto aa = static_cast<std::underlying_type_t<mutk::Potential::some_t>>(a);
        return pot.model_.CreateCountMatrix(n, t, aa);
    }

    template<class Arg>
    inline
    auto CreateGameteMatrix(size_t n, Arg arg, float t) const {
        auto ret = message_t::from_shape({
            message_axis_size(n, Ploidy::Diploid),
            message_axis_size(n, Ploidy::Haploid)
        });
        auto mat = CreateModelMatrix(n, arg, t);
        for(message_t::size_type i = 0; i < ret.shape(0); ++i) {
            for(message_t::size_type j = 0; j < ret.shape(1); ++j) {
                // a/b -> x
                auto [a,b] = mutk::diploid_alleles(i);
                auto x = mutk::haploid_allele(j);
                ret(i,j) = 0.5*(mat(a,x) + mat(b,x));
            }
        }
        return ret;
    }

    template<class Arg>
    inline
    auto CreateModelMatrixU(size_t n, Arg a) const {
        return CreateModelMatrix(n, a, pot.u_);
    }

    template<class Arg>
    inline
    auto CreateModelMatrixV(size_t n, Arg a) const {
        return CreateModelMatrix(n, a, pot.v_);
    }

    template<class Arg>
    inline
    auto CreateGameteMatrixU(size_t n, Arg a) const {
        return CreateGameteMatrix(n, a, pot.u_);
    }

    template<class Arg>
    inline
    auto CreateGameteMatrixV(size_t n, Arg a) const {
        return CreateGameteMatrix(n, a, pot.v_);
    }

    const mutk::MatingPotential &pot;
};

template<class Arg>
inline
message_t mutk::MatingPotential::Impl::operator()(size_t n, Arg a) {
    auto ploidy0 = message_axis_ploidy(pot.labels_.sequence()[0]);
    auto ploidy1 = message_axis_ploidy(pot.labels_.sequence()[1]);
    auto ploidy2 = message_axis_ploidy(pot.labels_.sequence()[2]);

    if(ploidy0 == Ploidy::Diploid) {
        if(ploidy1 == Ploidy::Diploid) {
            if(ploidy2 == Ploidy::Diploid) {
                return Create<222>::call(*this, n, a);
            } else {
                return Create<221>::call(*this, n, a);
            }
        } else if(ploidy2 == Ploidy::Diploid) {
            return Create<212>::call(*this, n, a);
        } else {
            return Create<211>::call(*this, n, a);
        }
    } else if(ploidy1 == Ploidy::Diploid) {
        if(ploidy2 == Ploidy::Diploid) {
            return Create<122>::call(*this, n, a);
        } else {
            return Create<121>::call(*this, n, a);
        }
    } else if(ploidy2 == Ploidy::Diploid) {
        return Create<112>::call(*this, n, a);
    } else {
        return Create<111>::call(*this, n, a);
    }
}

message_t mutk::MatingPotential::Create(size_t n, any_t a) {
    return MatingPotential::Impl(*this)(n,a);
}

message_t mutk::MatingPotential::Create(size_t n, mean_t a) {
    return MatingPotential::Impl(*this)(n,a);
}

message_t mutk::MatingPotential::Create(size_t n, some_t a) {
    return MatingPotential::Impl(*this)(n,a);
}

