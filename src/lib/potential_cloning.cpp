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

inline
auto create_model_matrix(const mutk::MutationModel &m, size_t n, float_t t, mutk::Potential::any_t) {
    return m.CreateTransitionMatrix(n, t);
}

inline
auto create_model_matrix(const mutk::MutationModel &m, size_t n, float_t t, mutk::Potential::mean_t a) {
    auto b = static_cast<std::underlying_type_t<mutk::Potential::mean_t>>(a);
    return (b > 0) ? m.CreateMeanMatrix(n, t) : m.CreateTransitionMatrix(n, t);
}

inline
auto create_model_matrix(const mutk::MutationModel &m, size_t n, float_t t, mutk::Potential::some_t a) {
    return m.CreateCountMatrix(n, t, +a);
}

template<>
template<class Arg>
message_t mutk::CloningPotential<Ploidy::Haploid, Ploidy::Haploid>::DoCreate(size_t n, Arg a) {
    return create_model_matrix(model_, u_, n, a);
}

template<>
template<class Arg>
message_t mutk::CloningPotential<Ploidy::Diploid, Ploidy::Haploid>::DoCreate(size_t n, Arg a) {
    auto ret = message_t::from_shape(MessageShape(n));
    auto mat = create_model_matrix(model_, u_, n, a);
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
message_t mutk::CloningPotential<Ploidy::Haploid, Ploidy::Diploid>::DoCreate(size_t n, Arg a) {
    auto ret = message_t::from_shape(MessageShape(n));
    auto mat = create_model_matrix(model_, u_, n, a);
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
//     - This works because `create_model_matrix(... mean_t(0))` returns the transition matrix.
template<>
template<class Arg>
message_t mutk::CloningPotential<Ploidy::Diploid, Ploidy::Diploid>::DoCreate(size_t n, Arg a) {
    using int_t = std::underlying_type_t<Arg>;
    auto ret = message_t::from_shape(MessageShape(n));
    ret.fill(0.0f);
    int_t count = static_cast<int_t>(a);
    for(int_t k=0; k<=count; ++k) {
        auto mat1 = create_model_matrix(model_, u_, n, Arg(k));
        auto mat2 = create_model_matrix(model_, u_, n, Arg(count-k));
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
TEST_CASE("CloningPotential<2,2>.Create") {
    using mutk::variable_t;

    using pot_t = mutk::CloningPotential<Ploidy::Diploid,Ploidy::Diploid>;

    mutk::Potential::labels_t labels({
        {variable_t{0}, Ploidy::Diploid},
        {variable_t{1}, Ploidy::Diploid}
    });

    auto test = [&](size_t n, float u, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = pot_t(labels, model, u);

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

        CHECK(matp.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(matpm.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(mat0.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(mat1.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(mat2.shape() == S{n*(n+1)/2,n*(n+1)/2});

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
TEST_CASE("CloningPotential<1,1>.Create") {
    using mutk::variable_t;

    using pot_t = mutk::CloningPotential<Ploidy::Haploid,Ploidy::Haploid>;

    mutk::Potential::labels_t labels({
        {variable_t{0}, Ploidy::Haploid},
        {variable_t{1}, Ploidy::Haploid}
    });

    auto test = [&](size_t n, float u, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = pot_t(labels, model, u);

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

        CHECK(matp.shape() == S{n,n});
        CHECK(matpm.shape() == S{n,n});
        CHECK(mat0.shape() == S{n,n});
        CHECK(mat1.shape() == S{n,n});
        CHECK(mat2.shape() == S{n,n});

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

// LCOV_EXCL_START
TEST_CASE("CloningPotential<2,1>.Create") {
    using mutk::variable_t;

    using pot_t = mutk::CloningPotential<Ploidy::Diploid,Ploidy::Haploid>;

    mutk::Potential::labels_t labels({
        {variable_t{0}, Ploidy::Diploid},
        {variable_t{1}, Ploidy::Haploid}
    });

    auto test = [&](size_t n, float u, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = pot_t(labels, model, u);

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

        CHECK(matp.shape() == S{n,n});
        CHECK(matpm.shape() == S{n,n});
        CHECK(mat0.shape() == S{n,n});
        CHECK(mat1.shape() == S{n,n});
        CHECK(mat2.shape() == S{n,n});

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
TEST_CASE("CloningPotential<1,2>.Create") {
    using mutk::variable_t;

    using pot_t = mutk::CloningPotential<Ploidy::Diploid,Ploidy::Haploid>;

    mutk::Potential::labels_t labels({
        {variable_t{0}, Ploidy::Diploid},
        {variable_t{1}, Ploidy::Haploid}
    });

    auto test = [&](size_t n, float u, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = pot_t(labels, model, u);

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

        CHECK(matp.shape() == S{n,n});
        CHECK(matpm.shape() == S{n,n});
        CHECK(mat0.shape() == S{n,n});
        CHECK(mat1.shape() == S{n,n});
        CHECK(mat2.shape() == S{n,n});

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
