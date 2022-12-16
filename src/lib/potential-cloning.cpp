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
struct mutk::CloningPotential::Impl {
    Impl(const mutk::CloningPotential &p) : pot{p} {}

    inline
    auto CreateModelMatrix(size_t n, mutk::Potential::any_t) const {
        return pot.model_.CreateTransitionMatrix(n, pot.u_);
    }

    inline
    auto CreateModelMatrix(size_t n, mutk::Potential::mean_t a) const {
        auto aa = static_cast<std::underlying_type_t<mutk::Potential::mean_t>>(a);
        return (aa > 0) ? pot.model_.CreateMeanMatrix(n, pot.u_) :
                          pot.model_.CreateTransitionMatrix(n, pot.u_);
    }

    inline
    auto CreateModelMatrix(size_t n, mutk::Potential::some_t a) const {
        auto aa = static_cast<std::underlying_type_t<mutk::Potential::some_t>>(a);
        return pot.model_.CreateCountMatrix(n, pot.u_, aa);
    }

    template<class Arg>
    inline
    message_t operator()(size_t n, Arg a);

    // Use nested class to mimic partial template specialization for functions
    template<int N>
    struct Create {
        template<class Arg>
        static message_t call(const Impl &impl, size_t n, Arg a);
    };

    const mutk::CloningPotential &pot;
};

template<class Arg>
inline
message_t mutk::CloningPotential::Impl::operator()(size_t n, Arg a) {
    auto ploidy0 = message_axis_ploidy(pot.labels_.sequence()[0]);
    auto ploidy1 = message_axis_ploidy(pot.labels_.sequence()[1]);

    if(ploidy0 == Ploidy::Diploid) {
        if(ploidy1 == Ploidy::Diploid) {
            return Create<22>::call(*this, n, a);
        } else {
            return Create<21>::call(*this, n, a);
        }
    } else if(ploidy1 == Ploidy::Diploid) {
        return Create<12>::call(*this, n, a);
    } else {
        return Create<11>::call(*this, n, a);
    }
}

message_t mutk::CloningPotential::Create(size_t n, any_t a) {
    return CloningPotential::Impl(*this)(n,a);
}

message_t mutk::CloningPotential::Create(size_t n, mean_t a) {
    return CloningPotential::Impl(*this)(n,a);
}

message_t mutk::CloningPotential::Create(size_t n, some_t a) {
    return CloningPotential::Impl(*this)(n,a);
}

template<>
template<class Arg>
message_t mutk::CloningPotential::Impl::Create<11>::call(const Impl &impl, size_t n, Arg a) {
    return impl.CreateModelMatrix(n, a);
}

template<>
template<class Arg>
message_t mutk::CloningPotential::Impl::Create<21>::call(const Impl &impl, size_t n, Arg a) {
    auto ret = message_t::from_shape(impl.pot.Shape(n));
    auto mat = impl.CreateModelMatrix(n, a);
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

template<>
template<class Arg>
message_t mutk::CloningPotential::Impl::Create<12>::call(const Impl &impl, size_t n, Arg a) {
    auto ret = message_t::from_shape(impl.pot.Shape(n));
    auto mat = impl.CreateModelMatrix(n, a);
    for(message_t::size_type i = 0; i < ret.shape(0); ++i) {
        for(message_t::size_type j = 0; j < ret.shape(1); ++j) {
            // a -> x/x
            auto a = mutk::haploid_allele(i);
            auto [x,y] = mutk::diploid_alleles(j);
            ret(i,j) = (x == y) ? mat(a,x) : 0.0f;
        }
    }
    return ret;
}

// If `a` is `any_t(0)` this will calculate P(j|i) via 1 loop.
// If `a` is `some_t(k)` this will calculate P(k mutations & j | i) via k+1 loops.
// If `a` is `mean_t(1)` this will calculate P(j|i)*E[number of mutations | j,i] via 2 loops.
//     - This works because `CreateModelMatrix(... mean_t(0))` returns the transition matrix.
template<>
template<class Arg>
message_t mutk::CloningPotential::Impl::Create<22>::call(const Impl &impl, size_t n, Arg a) {
    using int_t = std::underlying_type_t<Arg>;
    auto ret = message_t::from_shape(impl.pot.Shape(n));
    ret.fill(0.0f);
    int_t count = static_cast<int_t>(a);
    for(int_t k=0; k<=count; ++k) {
        auto mat1 = impl.CreateModelMatrix(n, Arg(k));
        auto mat2 = impl.CreateModelMatrix(n, Arg(count-k));
        for(message_size_t i = 0; i < ret.shape(0); ++i) {
            for(message_size_t j = 0; j < ret.shape(1); ++j) {
                // a/b -> x/y
                auto [a,b] = diploid_alleles(i);
                auto [x,y] = diploid_alleles(j);
                float_t value = mat1(a,x)*mat2(b,y);
                if(x != y) {
                    value += mat1(a,y)*mat2(b,x);
                }
                ret(i,j) += value;
            }
        }
    }
    return ret;
}

// Notes:
// Expectation Semiring:
// https://www.cs.jhu.edu/~jason/papers/li+eisner.emnlp09.pdf
//   <p,r>
//   Prod: <p1*p2, p1*r2 + p2*r1>
//   Sum: <p1 + p2, r1 + r2>
//   0: <0,0>
//   1: <1,0>

// LCOV_EXCL_START
TEST_CASE("CloningPotential.Create for Diploid-Diploid") {
    using mutk::variable_t;

    std::vector<mutk::message_label_t> labels({
        mutk::make_message_label(variable_t{0}, Ploidy::Diploid),
        mutk::make_message_label(variable_t{1}, Ploidy::Diploid)
    });

    auto test = [&](size_t n, float u, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = mutk::CloningPotential(model, u, labels);

        auto obs_any = pot.Create(n, mutk::Potential::ANY);
        auto obs_mean = pot.Create(n, mutk::Potential::MEAN);
        auto obs_zero = pot.Create(n, mutk::Potential::ZERO);
        auto obs_one = pot.Create(n, mutk::Potential::ONE);
        auto obs_two = pot.Create(n, mutk::Potential::TWO);

        auto matp = model.CreateTransitionMatrix(n, u);
        auto matpm = model.CreateMeanMatrix(n, u);
        auto mat0 = model.CreateCountMatrix(n, u, 0);
        auto mat1 = model.CreateCountMatrix(n, u, 1);
        auto mat2 = model.CreateCountMatrix(n, u, 2);

        CHECK(obs_any.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(obs_mean.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(obs_zero.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(obs_one.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(obs_two.shape() == S{n*(n+1)/2,n*(n+1)/2});

        for(size_t a1=0,a=0;a1<n;++a1) {
            for(size_t a2=0;a2<=a1;++a2,++a) {
                for(size_t b1=0,b=0;b1<n;++b1) {
                    for(size_t b2=0;b2<=b1;++b2,++b) {
                        // a1/a2 -> b1/b2
                        INFO("Transition: ", a1, "/", a2, " -> ", b1, "/", b2, " (", a, " -> ", b, ")");
                        // a1 -> b1 && a2 -> b2
                        float val_any  = matp(a1,b1)*matp(a2,b2);
                        float val_mean = matp(a1,b1)*matpm(a2,b2) + matpm(a1,b1)*matp(a2,b2);
                        float val_zero = mat0(a1,b1)*mat0(a2,b2);
                        float val_one =  mat1(a1,b1)*mat0(a2,b2) + mat0(a1,b1)*mat1(a2,b2);
                        float val_two =  mat2(a1,b1)*mat0(a2,b2) + mat1(a1,b1)*mat1(a2,b2) + mat0(a1,b1)*mat2(a2,b2);
                        if(b1 != b2) {
                            // a2 -> b1 && a1 -> b2
                            val_any  += matp(a2,b1)*matp(a1,b2);
                            val_mean += matp(a2,b1)*matpm(a1,b2) + matpm(a2,b1)*matp(a1,b2);
                            val_zero += mat0(a2,b1)*mat0(a1,b2);
                            val_one  += mat1(a2,b1)*mat0(a1,b2) + mat0(a2,b1)*mat1(a1,b2);
                            val_two  += mat2(a2,b1)*mat0(a1,b2) + mat1(a2,b1)*mat1(a1,b2) + mat0(a2,b1)*mat2(a1,b2);
                        }
                        CHECK(obs_any(a,b)  == doctest::Approx(val_any));
                        CHECK(obs_mean(a,b) == doctest::Approx(val_mean));
                        CHECK(obs_zero(a,b) == doctest::Approx(val_zero));
                        CHECK(obs_one(a,b)  == doctest::Approx(val_one));
                        CHECK(obs_two(a,b)  == doctest::Approx(val_two));
                    }
                }
            }
        }
    };
    run_mutation_tests(test);
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
TEST_CASE("CloningPotential.Create for Diploid-Haploid") {
    using mutk::variable_t;

    std::vector<mutk::message_label_t> labels({
        mutk::make_message_label(variable_t{0}, Ploidy::Diploid),
        mutk::make_message_label(variable_t{1}, Ploidy::Haploid)
    });
    auto test = [&](size_t n, float u, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = mutk::CloningPotential(model, u, labels);

        auto obs_any = pot.Create(n, mutk::Potential::ANY);
        auto obs_mean = pot.Create(n, mutk::Potential::MEAN);
        auto obs_zero = pot.Create(n, mutk::Potential::ZERO);
        auto obs_one = pot.Create(n, mutk::Potential::ONE);
        auto obs_two = pot.Create(n, mutk::Potential::TWO);

        auto matp = model.CreateTransitionMatrix(n, u);
        auto matpm = model.CreateMeanMatrix(n, u);
        auto mat0 = model.CreateCountMatrix(n, u, 0);
        auto mat1 = model.CreateCountMatrix(n, u, 1);
        auto mat2 = model.CreateCountMatrix(n, u, 2);

        CHECK(obs_any.shape() == S{n*(n+1)/2,n});
        CHECK(obs_mean.shape() == S{n*(n+1)/2,n});
        CHECK(obs_zero.shape() == S{n*(n+1)/2,n});
        CHECK(obs_one.shape() == S{n*(n+1)/2,n});
        CHECK(obs_two.shape() == S{n*(n+1)/2,n});

        for(size_t a1=0,a=0;a1<n;++a1) {
            for(size_t a2=0;a2<=a1;++a2,++a) {
                for(size_t b=0;b<n;++b) {
                    INFO("Transition: ", a1, "/", a2, " -> ", b, " (", a, " -> ", b, ")");
                    // a1/a2 -> b
                    float val_any  = 0.5f*(matp(a1,b)+matp(a2,b));
                    float val_mean = 0.5f*(matpm(a1,b)+matpm(a2,b));
                    float val_zero = 0.5f*(mat0(a1,b)+mat0(a2,b));
                    float val_one  = 0.5f*(mat1(a1,b)+mat1(a2,b));
                    float val_two  = 0.5f*(mat2(a1,b)+mat2(a2,b));

                    CHECK(obs_any(a,b)  == doctest::Approx(val_any));
                    CHECK(obs_mean(a,b) == doctest::Approx(val_mean));
                    CHECK(obs_zero(a,b) == doctest::Approx(val_zero));
                    CHECK(obs_one(a,b)  == doctest::Approx(val_one));
                    CHECK(obs_two(a,b)  == doctest::Approx(val_two));
                }
            }
        }
    };
    run_mutation_tests(test);
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
TEST_CASE("CloningPotential.Create for Haploid-Diploid") {
    using mutk::variable_t;

    std::vector<mutk::message_label_t> labels({
        mutk::make_message_label(variable_t{0}, Ploidy::Haploid),
        mutk::make_message_label(variable_t{1}, Ploidy::Diploid)
    });

    auto test = [&](size_t n, float u, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = mutk::CloningPotential(model, u, labels);

        auto obs_any = pot.Create(n, mutk::Potential::ANY);
        auto obs_mean = pot.Create(n, mutk::Potential::MEAN);
        auto obs_zero = pot.Create(n, mutk::Potential::ZERO);
        auto obs_one = pot.Create(n, mutk::Potential::ONE);
        auto obs_two = pot.Create(n, mutk::Potential::TWO);

        auto matp = model.CreateTransitionMatrix(n, u);
        auto matpm = model.CreateMeanMatrix(n, u);
        auto mat0 = model.CreateCountMatrix(n, u, 0);
        auto mat1 = model.CreateCountMatrix(n, u, 1);
        auto mat2 = model.CreateCountMatrix(n, u, 2);

        CHECK(obs_any.shape() == S{n,n*(n+1)/2});
        CHECK(obs_mean.shape() == S{n,n*(n+1)/2});
        CHECK(obs_zero.shape() == S{n,n*(n+1)/2});
        CHECK(obs_one.shape() == S{n,n*(n+1)/2});
        CHECK(obs_two.shape() == S{n,n*(n+1)/2});

        for(size_t a=0;a<n;++a) {
            for(size_t b1=0,b=0;b1<n;++b1) {
                for(size_t b2=0;b2<=b1;++b2,++b) {
                    INFO("Transition: ", a, " -> ", b1, "/", b2, " (", a, " -> ", b, ")");
                    // a -> b1/b2
                    float val_any  = (b1 == b2) ? matp(a,b1) : 0.0f;
                    float val_mean = (b1 == b2) ? matpm(a,b1) : 0.0f;
                    float val_zero = (b1 == b2) ? mat0(a,b1) : 0.0f;
                    float val_one  = (b1 == b2) ? mat1(a,b1) : 0.0f;
                    float val_two  = (b1 == b2) ? mat2(a,b1) : 0.0f;

                    CHECK(obs_any(a,b)  == doctest::Approx(val_any));
                    CHECK(obs_mean(a,b) == doctest::Approx(val_mean));
                    CHECK(obs_zero(a,b) == doctest::Approx(val_zero));
                    CHECK(obs_one(a,b)  == doctest::Approx(val_one));
                    CHECK(obs_two(a,b)  == doctest::Approx(val_two));
                }
            }
        }
    };
    run_mutation_tests(test);
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
TEST_CASE("CloningPotential.Create for Haploid-Haploid") {
    using mutk::variable_t;

    std::vector<mutk::message_label_t> labels({
        mutk::make_message_label(variable_t{0}, Ploidy::Haploid),
        mutk::make_message_label(variable_t{1}, Ploidy::Haploid)
    });

    auto test = [&](size_t n, float u, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = mutk::CloningPotential(model, u, labels);

        auto obs_any = pot.Create(n, mutk::Potential::ANY);
        auto obs_mean = pot.Create(n, mutk::Potential::MEAN);
        auto obs_zero = pot.Create(n, mutk::Potential::ZERO);
        auto obs_one = pot.Create(n, mutk::Potential::ONE);
        auto obs_two = pot.Create(n, mutk::Potential::TWO);

        auto matp = model.CreateTransitionMatrix(n, u);
        auto matpm = model.CreateMeanMatrix(n, u);
        auto mat0 = model.CreateCountMatrix(n, u, 0);
        auto mat1 = model.CreateCountMatrix(n, u, 1);
        auto mat2 = model.CreateCountMatrix(n, u, 2);

        CHECK(obs_any.shape() == S{n,n});
        CHECK(obs_mean.shape() == S{n,n});
        CHECK(obs_zero.shape() == S{n,n});
        CHECK(obs_one.shape() == S{n,n});
        CHECK(obs_two.shape() == S{n,n});

        for(size_t a=0;a<n;++a) {
            for(size_t b=0;b<n;++b) {
                INFO("Transition: ", a, " -> ", b, " (", a, " -> ", b, ")");
                // a -> b
                float val_any = matp(a,b);
                float val_mean = matpm(a,b);
                float val_zero = mat0(a,b);
                float val_one = mat1(a,b);
                float val_two = mat2(a,b);

                CHECK(obs_any(a,b)  == doctest::Approx(val_any));
                CHECK(obs_mean(a,b) == doctest::Approx(val_mean));
                CHECK(obs_zero(a,b) == doctest::Approx(val_zero));
                CHECK(obs_one(a,b)  == doctest::Approx(val_one));
                CHECK(obs_two(a,b)  == doctest::Approx(val_two));
            }
        }
    };
    run_mutation_tests(test);
}
// LCOV_EXCL_STOP
