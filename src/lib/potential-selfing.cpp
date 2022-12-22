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
struct mutk::SelfingPotential::Impl {
    Impl(const mutk::SelfingPotential &p) : pot{p} {}

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

    template<class Arg>
    inline
    message_t operator()(size_t n, Arg a);

    // Use nested class to mimic partial template specialization for functions
    template<int N>
    struct Create {
        template<class Arg>
        static message_t call(const Impl &impl, size_t n, Arg a);
    };

    const mutk::SelfingPotential &pot;
};

template<class Arg>
inline
message_t mutk::SelfingPotential::Impl::operator()(size_t n, Arg a) {
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

message_t mutk::SelfingPotential::Create(size_t n, any_t a) {
    return SelfingPotential::Impl(*this)(n,a);
}

message_t mutk::SelfingPotential::Create(size_t n, mean_t a) {
    return SelfingPotential::Impl(*this)(n,a);
}

message_t mutk::SelfingPotential::Create(size_t n, some_t a) {
    return SelfingPotential::Impl(*this)(n,a);
}

template<>
template<class Arg>
message_t mutk::SelfingPotential::Impl::Create<22>::call(const Impl &impl, size_t n, Arg arg) {
    using int_t = std::underlying_type_t<Arg>;
    auto ret = message_t::from_shape(impl.pot.Shape(n));
    ret.fill(0.0f);
    int_t val = static_cast<int_t>(arg);
    for(int_t k=0; k<=val; ++k) {
        auto mat1 = impl.CreateModelMatrixU(n, Arg(k));
        auto mat2 = impl.CreateModelMatrixV(n, Arg(val-k));
        for(message_size_t i = 0; i < ret.shape(0); ++i) {
            for(message_size_t j = 0; j < ret.shape(1); ++j) {
                // (a/b -> x) * (a/b -> y)  -> x/y
                auto [a,b] = diploid_alleles(i);
                auto [x,y] = diploid_alleles(j);
                float_t value = (mat1(a,x)+mat1(b,x))*(mat2(a,y)+mat2(b,y));
                if(x != y) {
                    value += (mat1(a,y)+mat1(b,y))*(mat2(a,x)+mat2(b,x));
                }
                ret(i,j) += 0.25*value;
            }
        }
    }
    return ret;
}

template<>
template<class Arg>
message_t mutk::SelfingPotential::Impl::Create<12>::call(const Impl &impl, size_t n, Arg arg) {
    // A haploid individual selfing (1x1 -> 2) is equivalent to generating two clones with
    // different mutation rates and those clones fusing.
    using int_t = std::underlying_type_t<Arg>;
    auto ret = message_t::from_shape(impl.pot.Shape(n));
    ret.fill(0.0f);
    int_t val = static_cast<int_t>(arg);
    for(int_t k=0; k<=val; ++k) {
        auto mat1 = impl.CreateModelMatrixU(n, Arg(k));
        auto mat2 = impl.CreateModelMatrixV(n, Arg(val-k));
        for(message_size_t i = 0; i < ret.shape(0); ++i) {
            for(message_size_t j = 0; j < ret.shape(1); ++j) {
                // a -> x/y
                auto a = haploid_allele(i);
                auto [x,y] = diploid_alleles(j);
                float_t value = mat1(a,x)*mat2(a,y);
                if(x != y) {
                    value += mat1(a,y)*mat2(a,x);
                }
                ret(i,j) += value;
            }
        }
    }
    return ret;
}

template<>
template<class Arg>
message_t mutk::SelfingPotential::Impl::Create<21>::call(const Impl &impl, size_t n, Arg arg) {
    auto ret = message_t::from_shape(impl.pot.Shape(n));
    auto matU = impl.CreateModelMatrixU(n, arg);
    auto matV = impl.CreateModelMatrixV(n, arg);
    for(message_t::size_type i = 0; i < ret.shape(0); ++i) {
        for(message_t::size_type j = 0; j < ret.shape(1); ++j) {
            // a/b -> x
            auto [a,b] = mutk::diploid_alleles(i);
            auto x = mutk::haploid_allele(j);
            ret(i,j) = 0.25*(matU(a,x) + matU(b,x) + matV(a,x) + matV(b,x));
        }
    }
    return ret;
}

template<>
template<class Arg>
message_t mutk::SelfingPotential::Impl::Create<11>::call(const Impl &impl, size_t n, Arg arg) {
    // A haploid individual selfing and producing a haploid offspring is equivalent to:
    //   (1x1 -> 2 -> 1)
    //   The haploid individual generates two clones with different mutation rates.
    //   The clone fuse. Then the fused individual produces a haploid gamete.
    auto matU = impl.CreateModelMatrixU(n, arg);
    auto matV = impl.CreateModelMatrixV(n, arg);
    return 0.5*(matU+matV);
}

// LCOV_EXCL_START
TEST_CASE("SelfingPotential.Create for Diploid-Diploid") {
    using mutk::variable_t;

    std::vector<mutk::message_label_t> labels({
        mutk::make_message_label(variable_t{0}, Ploidy::Diploid),
        mutk::make_message_label(variable_t{1}, Ploidy::Diploid)
    });

    auto test = [&](size_t n, float u, float v, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u); CAPTURE(v);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = mutk::SelfingPotential(model, u, v, labels);

        auto obs_any = pot.Create(n, mutk::Potential::ANY);
        auto obs_mean = pot.Create(n, mutk::Potential::MEAN);
        auto obs_zero = pot.Create(n, mutk::Potential::ZERO);
        auto obs_one = pot.Create(n, mutk::Potential::ONE);
        auto obs_two = pot.Create(n, mutk::Potential::TWO);

        auto umatp = model.CreateTransitionMatrix(n, u);
        auto umatpm = model.CreateMeanMatrix(n, u);
        auto umat0 = model.CreateCountMatrix(n, u, 0);
        auto umat1 = model.CreateCountMatrix(n, u, 1);
        auto umat2 = model.CreateCountMatrix(n, u, 2);

        auto vmatp = model.CreateTransitionMatrix(n, v);
        auto vmatpm = model.CreateMeanMatrix(n, v);
        auto vmat0 = model.CreateCountMatrix(n, v, 0);
        auto vmat1 = model.CreateCountMatrix(n, v, 1);
        auto vmat2 = model.CreateCountMatrix(n, v, 2);

        CHECK(obs_any.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(obs_mean.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(obs_zero.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(obs_one.shape() == S{n*(n+1)/2,n*(n+1)/2});
        CHECK(obs_two.shape() == S{n*(n+1)/2,n*(n+1)/2});

        for(size_t a1=0,a=0;a1<n;++a1) {
            for(size_t a2=0;a2<=a1;++a2,++a) {
                size_t b1 = a1;
                size_t b2 = a2;
                for(size_t c1=0,c=0;c1<n;++c1) {
                    for(size_t c2=0;c2<=c1;++c2,++c) {
                        // a1/a2 -> c1/c2 as a1/a2 -> c1 && b1/b2 -> c2
                        INFO("Transition: ", a1, "/", a2, "*", b1, "/", b2, " -> ", c1, "/", c2, " (", a, " -> ", c, ")");

                        float val_any = 0.0;
                        float val_mean = 0.0;
                        float val_zero = 0.0;
                        float val_one = 0.0;
                        float val_two = 0.0;
                        for(auto aa : {a1,a2}) {
                            for(auto bb : {b1,b2}) {
                                // aa -> c1 && bb -> c2
                                val_any  += 0.25*umatp(aa,c1)*vmatp(bb,c2);
                                val_mean += 0.25*umatp(aa,c1)*vmatpm(bb,c2) + 0.25*umatpm(aa,c1)*vmatp(bb,c2);
                                val_zero += 0.25*umat0(aa,c1)*vmat0(bb,c2);
                                val_one  += 0.25*umat1(aa,c1)*vmat0(bb,c2)+0.25*umat0(aa,c1)*vmat1(bb,c2);
                                val_two  += 0.25*umat2(aa,c1)*vmat0(bb,c2)+0.25*umat1(aa,c1)*vmat1(bb,c2)+
                                            0.25*umat0(aa,c1)*vmat2(bb,c2);
                                if(c1 != c2) {
                                    // aa -> c2 && bb -> c1
                                    val_any  += 0.25*umatp(aa,c2)*vmatp(bb,c1);
                                    val_mean += 0.25*umatp(aa,c2)*vmatpm(bb,c1) + 0.25*umatpm(aa,c2)*vmatp(bb,c1);
                                    val_zero += 0.25*umat0(aa,c2)*vmat0(bb,c1);
                                    val_one  += 0.25*umat1(aa,c2)*vmat0(bb,c1)+0.25*umat0(aa,c2)*vmat1(bb,c1);
                                    val_two  += 0.25*umat2(aa,c2)*vmat0(bb,c1)+0.25*umat1(aa,c2)*vmat1(bb,c1)+
                                                0.25*umat0(aa,c2)*vmat2(bb,c1);
                                }
                            }
                        }
                        CHECK(obs_any(a,c)  == doctest::Approx(val_any));
                        CHECK(obs_mean(a,c) == doctest::Approx(val_mean));
                        CHECK(obs_zero(a,c) == doctest::Approx(val_zero));
                        CHECK(obs_one(a,c)  == doctest::Approx(val_one));
                        CHECK(obs_two(a,c)  == doctest::Approx(val_two));
                    }
                }
            }
        }
    };
    run_mutation_tests(test);
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
TEST_CASE("SelfingPotential.Create for Diploid-Haploid") {
    using mutk::variable_t;

    std::vector<mutk::message_label_t> labels({
        mutk::make_message_label(variable_t{0}, Ploidy::Diploid),
        mutk::make_message_label(variable_t{1}, Ploidy::Haploid)
    });

    auto test = [&](size_t n, float u, float v, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u); CAPTURE(v);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = mutk::SelfingPotential(model, u, v, labels);

        auto obs_any = pot.Create(n, mutk::Potential::ANY);
        auto obs_mean = pot.Create(n, mutk::Potential::MEAN);
        auto obs_zero = pot.Create(n, mutk::Potential::ZERO);
        auto obs_one = pot.Create(n, mutk::Potential::ONE);
        auto obs_two = pot.Create(n, mutk::Potential::TWO);

        auto umatp = model.CreateTransitionMatrix(n, u);
        auto umatpm = model.CreateMeanMatrix(n, u);
        auto umat0 = model.CreateCountMatrix(n, u, 0);
        auto umat1 = model.CreateCountMatrix(n, u, 1);
        auto umat2 = model.CreateCountMatrix(n, u, 2);

        auto vmatp = model.CreateTransitionMatrix(n, v);
        auto vmatpm = model.CreateMeanMatrix(n, v);
        auto vmat0 = model.CreateCountMatrix(n, v, 0);
        auto vmat1 = model.CreateCountMatrix(n, v, 1);
        auto vmat2 = model.CreateCountMatrix(n, v, 2);

        CHECK(obs_any.shape() == S{n*(n+1)/2,n});
        CHECK(obs_mean.shape() == S{n*(n+1)/2,n});
        CHECK(obs_zero.shape() == S{n*(n+1)/2,n});
        CHECK(obs_one.shape() == S{n*(n+1)/2,n});
        CHECK(obs_two.shape() == S{n*(n+1)/2,n});

        for(size_t a1=0,a=0;a1<n;++a1) {
            for(size_t a2=0;a2<=a1;++a2,++a) {
                size_t b1 = a1;
                size_t b2 = a2;
                for(size_t c=0; c < n; ++c) {
                    // a1/a2 -> c as a1/a2 -> c1 && b1/b2 -> c2 && c1 or c2 -> c
                    INFO("Transition: ", a1, "/", a2, "*", b1, "/", b2, " -> ", c, " (", a, " -> ", c, ")");
                    double val_any = 0.0;
                    double val_mean = 0.0;
                    double val_zero = 0.0;
                    double val_one = 0.0;
                    double val_two = 0.0;
                    for(auto aa : {a1,a2}) {
                        for(auto bb : {b1,b2}) {
                            // aa -> c || bb -> c
                            val_any  += 0.25*0.5*(umatp(aa,c)+vmatp(bb,c));
                            val_mean += 0.25*0.5*(umatpm(aa,c)+vmatpm(bb,c));
                            val_zero += 0.25*0.5*(umat0(aa,c)+vmat0(bb,c));
                            val_one  += 0.25*0.5*(umat1(aa,c)+vmat1(bb,c));
                            val_two  += 0.25*0.5*(umat2(aa,c)+vmat2(bb,c));
                        }
                    }
                    CHECK(obs_any(a,c)  == doctest::Approx(val_any));
                    CHECK(obs_mean(a,c) == doctest::Approx(val_mean));
                    CHECK(obs_zero(a,c) == doctest::Approx(val_zero));
                    CHECK(obs_one(a,c)  == doctest::Approx(val_one));
                    CHECK(obs_two(a,c)  == doctest::Approx(val_two));
                }
            }
        }
    };
    run_mutation_tests(test);
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
TEST_CASE("SelfingPotential.Create for Haploid-Diploid") {
    using mutk::variable_t;

    std::vector<mutk::message_label_t> labels({
        mutk::make_message_label(variable_t{0}, Ploidy::Haploid),
        mutk::make_message_label(variable_t{1}, Ploidy::Diploid)
    });

    auto test = [&](size_t n, float u, float v, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u); CAPTURE(v);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = mutk::SelfingPotential(model, u, v, labels);

        auto obs_any = pot.Create(n, mutk::Potential::ANY);
        auto obs_mean = pot.Create(n, mutk::Potential::MEAN);
        auto obs_zero = pot.Create(n, mutk::Potential::ZERO);
        auto obs_one = pot.Create(n, mutk::Potential::ONE);
        auto obs_two = pot.Create(n, mutk::Potential::TWO);

        auto umatp = model.CreateTransitionMatrix(n, u);
        auto umatpm = model.CreateMeanMatrix(n, u);
        auto umat0 = model.CreateCountMatrix(n, u, 0);
        auto umat1 = model.CreateCountMatrix(n, u, 1);
        auto umat2 = model.CreateCountMatrix(n, u, 2);

        auto vmatp = model.CreateTransitionMatrix(n, v);
        auto vmatpm = model.CreateMeanMatrix(n, v);
        auto vmat0 = model.CreateCountMatrix(n, v, 0);
        auto vmat1 = model.CreateCountMatrix(n, v, 1);
        auto vmat2 = model.CreateCountMatrix(n, v, 2);

        CHECK(obs_any.shape() == S{n,n*(n+1)/2});
        CHECK(obs_mean.shape() == S{n,n*(n+1)/2});
        CHECK(obs_zero.shape() == S{n,n*(n+1)/2});
        CHECK(obs_one.shape() == S{n,n*(n+1)/2});
        CHECK(obs_two.shape() == S{n,n*(n+1)/2});

        for(size_t a=0;a<n;++a) {
            size_t b = a;
            for(size_t c1=0,c=0;c1<n;++c1) {
                for(size_t c2=0;c2<=c1;++c2,++c) {
                    // a -> c1/c2 as a -> c1 && b -> c2
                    INFO("Transition: ", a, "*", b, " -> ", c1, "/", c2, " (", a, " -> ", c, ")");
                    double val_any = 0.0;
                    double val_mean = 0.0;
                    double val_zero = 0.0;
                    double val_one = 0.0;
                    double val_two = 0.0;
                    for(auto aa : {a}) {
                        for(auto bb : {b}) {
                            // aa -> c1 && bb -> c2
                            val_any  += umatp(aa,c1)*vmatp(bb,c2);
                            val_mean += umatp(aa,c1)*vmatpm(bb,c2) + umatpm(aa,c1)*vmatp(bb,c2);
                            val_zero += umat0(aa,c1)*vmat0(bb,c2);
                            val_one  += umat1(aa,c1)*vmat0(bb,c2) + umat0(aa,c1)*vmat1(bb,c2);
                            val_two  += umat2(aa,c1)*vmat0(bb,c2) + umat1(aa,c1)*vmat1(bb,c2)+
                                        umat0(aa,c1)*vmat2(bb,c2);
                            if(c1 != c2) {
                                // aa -> c2 && bb -> c1
                                val_any  += umatp(aa,c2)*vmatp(bb,c1);
                                val_mean += umatp(aa,c2)*vmatpm(bb,c1) + umatpm(aa,c2)*vmatp(bb,c1);
                                val_zero += umat0(aa,c2)*vmat0(bb,c1);
                                val_one  += umat1(aa,c2)*vmat0(bb,c1) + umat0(aa,c2)*vmat1(bb,c1);
                                val_two  += umat2(aa,c2)*vmat0(bb,c1) + umat1(aa,c2)*vmat1(bb,c1)+
                                            umat0(aa,c2)*vmat2(bb,c1);
                            }
                        }
                    }
                    CHECK(obs_any(a,c)  == doctest::Approx(val_any));
                    CHECK(obs_mean(a,c) == doctest::Approx(val_mean));
                    CHECK(obs_zero(a,c) == doctest::Approx(val_zero));
                    CHECK(obs_one(a,c)  == doctest::Approx(val_one));
                    CHECK(obs_two(a,c)  == doctest::Approx(val_two));
                }
            }
        }
    };
    run_mutation_tests(test);
}
// LCOV_EXCL_STOP

TEST_CASE("SelfingPotential.Create for Haploid-Haploid") {
    using mutk::variable_t;

    std::vector<mutk::message_label_t> labels({
        mutk::make_message_label(variable_t{0}, Ploidy::Haploid),
        mutk::make_message_label(variable_t{1}, Ploidy::Haploid)
    });

    auto test = [&](size_t n, float u, float v, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u); CAPTURE(v);

        using S = mutk::message_shape_t;
        S{0,0}; // silences a warning

        mutk::MutationModel model(k, 0.001, 0, 0, 0);
        auto pot = mutk::SelfingPotential(model, u, v, labels);

        auto obs_any = pot.Create(n, mutk::Potential::ANY);
        auto obs_mean = pot.Create(n, mutk::Potential::MEAN);
        auto obs_zero = pot.Create(n, mutk::Potential::ZERO);
        auto obs_one = pot.Create(n, mutk::Potential::ONE);
        auto obs_two = pot.Create(n, mutk::Potential::TWO);

        auto umatp = model.CreateTransitionMatrix(n, u);
        auto umatpm = model.CreateMeanMatrix(n, u);
        auto umat0 = model.CreateCountMatrix(n, u, 0);
        auto umat1 = model.CreateCountMatrix(n, u, 1);
        auto umat2 = model.CreateCountMatrix(n, u, 2);

        auto vmatp = model.CreateTransitionMatrix(n, v);
        auto vmatpm = model.CreateMeanMatrix(n, v);
        auto vmat0 = model.CreateCountMatrix(n, v, 0);
        auto vmat1 = model.CreateCountMatrix(n, v, 1);
        auto vmat2 = model.CreateCountMatrix(n, v, 2);

        CHECK(obs_any.shape() == S{n,n});
        CHECK(obs_mean.shape() == S{n,n});
        CHECK(obs_zero.shape() == S{n,n});
        CHECK(obs_one.shape() == S{n,n});
        CHECK(obs_two.shape() == S{n,n});

        for(size_t a=0;a<n;++a) {
            size_t b = a;
            for(size_t c=0;c<n;++c) {
                // a -> c as a -> c1 && b -> c2 && c1/c2 -> c
                INFO("Transition: ", a, "*", b, " -> ", c, " (", a, " -> ", c, ")");
                double val_any = 0.0;
                double val_mean = 0.0;
                double val_zero = 0.0;
                double val_one = 0.0;
                double val_two = 0.0;
                for(auto aa : {a}) {
                    for(auto bb : {b}) {
                        // aa -> c1 && bb -> c2
                        val_any  += 0.5*(umatp(aa,c)+vmatp(bb,c));
                        val_mean += 0.5*(umatpm(aa,c)+vmatpm(bb,c));
                        val_zero += 0.5*(umat0(aa,c)+vmat0(bb,c));
                        val_one  += 0.5*(umat1(aa,c)+vmat1(bb,c));
                        val_two  += 0.5*(umat2(aa,c)+vmat2(bb,c));
                    }
                }
                CHECK(obs_any(a,c)  == doctest::Approx(val_any));
                CHECK(obs_mean(a,c) == doctest::Approx(val_mean));
                CHECK(obs_zero(a,c) == doctest::Approx(val_zero));
                CHECK(obs_one(a,c)  == doctest::Approx(val_one));
                CHECK(obs_two(a,c)  == doctest::Approx(val_two));
            }
        }
    };
    run_mutation_tests(test);
}
// LCOV_EXCL_STOP
