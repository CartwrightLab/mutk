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
#include "unit_testing.hpp"
#include "mutation_testing.hpp"

#include <mutk/mutation.hpp>
#include <mutk/relationship_graph.hpp>

namespace {
constexpr int ALLELE[][2] = {
    {0,0},
    {0,1}, {1,1},
    {0,2}, {1,2}, {2,2},
    {0,3}, {1,3}, {2,3}, {3,3},
    {0,4}, {1,4}, {2,4}, {3,4}, {4,4},
    {0,5}, {1,5}, {2,5}, {3,5}, {4,5}, {5,5}
};
}

using mutk::MutationModel;
using mutk::message_t;

// ret(i,j) = P(j|i)
MutationModel::array_t MutationModel::CreateTransitionMatrix(message_size_t n, float_t t) const {
    assert(n > 0);
    assert(n <= 5);

    double beta = k_/(k_-1.0);
    double p_ij = -1.0/k_*expm1(-beta*t);
    double p_ii = exp(-beta*t) + p_ij;

    array_t ret = array_t::from_shape({n,n});

    for(message_size_t i = 0; i < n; ++i) {
        for(message_size_t j = 0; j < n; ++j) {
            ret(i,j) = (i == j) ? p_ii : p_ij;
        }
    }
    return ret;
}

// ret(i,j) = P(j & x mutations | i)
//
// beta = k/(k-1)
// P(x mutations | i) = (t^x Exp[-t])/x!
// 
MutationModel::array_t MutationModel::CreateCountMatrix(message_size_t n, float_t t, int x) const {
    assert(n > 0);
    assert(n <= 5);
    assert(x >= 0);
    double xlogt = (x == 0 && t == 0.0) ? 0.0 : x*log(t);
    double p_x = exp(-t+xlogt-lgamma(x+1));
    
    // Formula calculated using MatrixPower[] in Mathematica
    // Can be proved using an induction proof.
    double h = k_-1.0;
    double a = pow(-1.0/h, x);
    double p_ii = (1.0+h*a)/k_;
    double p_ij = (1.0-a)/k_;

    array_t ret = array_t::from_shape({n,n});

    for(size_t i = 0; i < n; ++i) {
        for(size_t j = 0; j < n; ++j) {
            ret(i,j) = (i == j) ? p_x*p_ii : p_x*p_ij;
        }
    }
    return ret;
}

// ret(i,j) = E[num of mutations | i,j]*P(j|i)
// Estimated via Mathematica
MutationModel::array_t MutationModel::CreateMeanMatrix(message_size_t n, float_t t) const {
    assert(n > 0);
    assert(n <= 5);

    double beta = k_/(k_-1.0);
    double p = -expm1(-beta*t);

    double p_ii = p * t/k_;
    double p_ij = (k_-p)/(k_-1.0) * t/k_;

    array_t ret = array_t::from_shape({n,n});

    for(message_size_t i = 0; i < n; ++i) {
        for(message_size_t j = 0; j < n; ++j) {
            ret(i,j) = (i == j) ? p_ii : p_ij;
        }
    }
    return ret;
}

// LCOV_EXCL_START
TEST_CASE("MutationModel.Constructor") {
    CHECK_NOTHROW(MutationModel(4.0, 0.001, 0.0, 0.0, 0.0));
    CHECK_THROWS_AS(MutationModel(1.0, 0.001, 0.0, 0.0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(MutationModel(4.0, -0.1, 0.0, 0.0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(MutationModel(4.0, 0.001, 1.1, 0.0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(MutationModel(4.0, 0.001, -3000, 0.0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(MutationModel(4.0, 0.001, 0.0, 1.1, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(MutationModel(4.0, 0.001, 0.0, -4000, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(MutationModel(4.0, 0.001, 0.0, 0.0, 1.1), std::invalid_argument);
    CHECK_THROWS_AS(MutationModel(4.0, 0.001, 0.0, 0.0, -3000), std::invalid_argument);
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
TEST_CASE("MutationModel.CreateTransitionMatrix") {
    using namespace boost::numeric::ublas;
    
    auto test = [&](size_t n, float u, float, float k) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

       MutationModel model(k, 0.001, 0, 0, 0);

        using mat_t = kalleles_test_mat::mat_t;

        kalleles_test_mat mat(n,k);

        mat_t P = identity_matrix<float>(n+1);
        mat_t m = P;
        float t = 1.0;
        // Use the first 11 elements of the expm series 
        for(int g=1;g <= 11; ++g) {
            mat_t a = prec_prod(m,mat.Q);
            m = a;
            t *= (u/g);
            P += m * t;
        }
        auto obs = model.CreateTransitionMatrix(n, u);
        REQUIRE(obs.dimension() == 2);
        REQUIRE(obs.shape(0) == n);
        REQUIRE(obs.shape(1) == n);
        for(size_t i=0; i < n; ++i) {
            for(size_t j=0; j < n; ++j) {
                INFO("Transition: ", i, " -> ", j);
                CHECK(obs(i,j) == doctest::Approx(P(i,j)));
            }
        }
    };
    run_mutation_tests(test);
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
TEST_CASE("MutationModel.CreateMeanMatrix") {
    using namespace boost::numeric::ublas;
    
    // TODO fix this test
    auto test = [&](size_t n, float u, float, float k) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        MutationModel model(k, 0.001, 0, 0, 0);

        using mat_t = kalleles_test_mat::mat_t;

        kalleles_test_mat mat(n,k);

        // First 11 elements of Poisson:
        // u^x*exp(-u)/x!
        // m = (Q+I)^x
        mat_t m = identity_matrix<float>(n+1);
        mat_t J = mat.Q+m;
        float t = exp(-u);
        mat_t P = m*t;
        mat_t S = P*0.0;
        for(int x=1; x <= 10; ++x) {
            mat_t a = prec_prod(J,m);
            m = a;
            t *= u/x;
            P = m*t;
            S += P*x;
        }
        auto obs = model.CreateMeanMatrix(n, u);
        REQUIRE(obs.dimension() == 2);
        REQUIRE(obs.shape(0) == n);
        REQUIRE(obs.shape(1) == n);
        for(size_t i=0; i < n; ++i) {
            for(size_t j=0; j < n; ++j) {
                INFO("Transition: ", i, " -> ", j);
                CHECK(obs(i,j) == doctest::Approx(S(i,j)));
            }
        }
    };

    run_mutation_tests(test);
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
TEST_CASE("MutationModel.CreateCountMatrix") {
    using namespace boost::numeric::ublas;
    
    auto test = [&](size_t n, float u, float, float k) {
        CAPTURE(n); CAPTURE(k); CAPTURE(u);

        MutationModel model(k, 0.001, 0, 0, 0);

        using mat_t = kalleles_test_mat::mat_t;

        kalleles_test_mat mat(n,k);

        mat_t P = identity_matrix<float>(n+1);
        mat_t m = P;
        mat_t J = mat.Q+m;
        float t;
        for(int x=0; x <= 10; ++x) {
            CAPTURE(x);
            if(x == 0) {
                P = m*exp(-u);
                t = 1.0;
            } else {
                mat_t a = prec_prod(J,m);
                m = a;
                t *= (u/x);
                P = m*(exp(-u)*t);
            }
            auto obs = model.CreateCountMatrix(n, u, x);
            REQUIRE(obs.dimension() == 2);
            REQUIRE(obs.shape(0) == n);
            REQUIRE(obs.shape(1) == n);
            for(size_t i=0; i < n; ++i) {
                for(size_t j=0; j < n; ++j) {
                    INFO("Transition: ", i, " -> ", j);
                    CAPTURE(i);
                    CAPTURE(j);
                    CHECK(obs(i,j) == doctest::Approx(P(i,j)));
                }
            }

        }
    };

    run_mutation_tests(test);
}
// LCOV_EXCL_STOP

#if 0

KAllelesModel::tensor_t KAllelesModel::CreatePriorHaploid(size_t n) const {
    double k = k_;
    double e = theta_/(k-1.0);

    double p_R = (1.0+e+(k-1.0)*e*hap_bias_)/(1.0+k*e);
    double p_A = (e-e*hap_bias_)/(1.0+k*e);

    tensor_t ret = tensor_t::from_shape({n});

    for(size_t i = 0; i < n; ++i) {
        ret(i) = (i == 0) ? p_R : p_A;
    }
    return ret;
}

// LCOV_EXCL_START
TEST_CASE("KAllelesModel-CreatePriorHaploid") {
    auto test_haploid = [](size_t n, float theta, float hap_bias,
        float k) {
        CAPTURE(n);
        CAPTURE(theta);
        CAPTURE(hap_bias);
        CAPTURE(k);

        KAllelesModel model(k, theta, 0, 0, hap_bias);

        auto obs = model.CreatePriorHaploid(n);

        if(n == k) {
            float s = xt::sum(obs)();
            CHECK(s == doctest::Approx(1.0f));
        }

        const int sz = n;
        for(int a=0;a<sz;++a) {
            CAPTURE(a);
            float e = theta/(k-1.0);
            float r;
            if(a == 0) {
                // reference haploid
                r = (1.0+e+e*(k-1)*hap_bias)/(1.0+k*e);
            } else {
                r = (e-e*hap_bias)/(1.0+k*e);
            }
            float expected = r;
            CHECK(obs(a) == doctest::Approx(expected));
        }
    };

    auto test = [&](float theta, float hap_bias, float k) {
        for(int i=1;i<=5;++i) {
            test_haploid(i, theta, hap_bias, k);
        }
    };

    test(0.001, 0.0, 5.0);
    test(0.01, 1.0, 6.0);
    test(0.1, 0.0, 4.0);
    test(0.001, 1.0, 5.0);
    test(0.001, -1.0, 5.5);
}
// LCOV_EXCL_STOP

KAllelesModel::tensor_t KAllelesModel::CreatePriorDiploid(size_t n) const {
    double k = k_;
    double e = theta_/(k-1.0);

    double p_hom = (1.0+e)/(1.0+k*e);
    double p_hetk = e/(1.0+k*e);

    double p_RR = p_hom*(2.0+e+(k-1.0)*e*hom_bias_)/(2.0+k*e);
    double p_AA = p_hom*(e-e*hom_bias_)/(2.0+k*e);

    double p_RA = p_hetk*(2.0+2.0*e+(k-2.0)*e*het_bias_)/(2.0+k*e);
    double p_AB = p_hetk*(2.0*e-2.0*e*het_bias_)/(2.0+k*e);

    tensor_t ret = tensor_t::from_shape({mutk::dim_width<2>(n)});

    for(size_t i = 0; i < ret.size(); ++i) {
        auto a = ALLELE[i][0];
        auto b = ALLELE[i][1];
        if(a == b) {
            ret(i) = (a == 0) ? p_RR : p_AA;
        } else {
            ret(i) = (a == 0) ? p_RA : p_AB;
        }
    }
    return ret;
}

// LCOV_EXCL_START
TEST_CASE("KAllelesModel-CreatePriorDiploid") {
    auto test_diploid = [](size_t n, float theta, float hom_bias,
        float het_bias, float k) {
        CAPTURE(n);
        CAPTURE(theta);
        CAPTURE(hom_bias);
        CAPTURE(het_bias);
        CAPTURE(k);

        KAllelesModel model(k, theta, hom_bias, het_bias, 0);

        auto obs = model.CreatePriorDiploid(n);

        if(n == k) {
            float s = xt::sum(obs)();
            CHECK(s == doctest::Approx(1.0f));
        }

        const int sz = n;
        int y = 0;
        for(int a=0;a<sz;++a) {
            for(int b=0;b<=a;++b) {
                CAPTURE(a);
                CAPTURE(b);
                float e = theta/(k-1.0);

                double r = 0.0;
                if(a == 0 && b == 0) {
                    // Reference Homozygote
                    r = (1.0+e)/(1.0+k*e)*(2.0+e+(k-1.0)*e*hom_bias)/(2.0+k*e);
                } else if(a == b) {
                    r = (1.0+e)/(1.0+k*e)*(e-e*hom_bias)/(2.0+k*e);
                } else if(b == 0 || a == 0) {
                    // Reference Heterozygote
                    r = e/(1.0+k*e)*(2.0+2.0*e+(k-2.0)*e*het_bias)/(2.0+k*e);
                } else {
                    // Alt heterozygote
                    r = e/(1.0+k*e)*(2.0*e-2.0*e*het_bias)/(2.0+k*e);
                }
                float expected = r;
                CHECK(obs(y) == doctest::Approx(expected));
                ++y;
            }
        }
    };

    auto test = [&](float theta, float hom_bias, float het_bias, float k) {
        for(int i=1;i<=5;++i) {
            test_diploid(i, theta, hom_bias, het_bias, k);
        }
    };

    test(0.001, 0.0, 0.0, 5.0);
    test(0.01, 1.0, 0.0, 6.0);
    test(0.1, 0.0, 1.0, 4.0);
    test(0.001, 1.0, 1.0, 5.0);
    test(0.001, -1.0, -1.0, 5.5);
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_clone_haploid_impl(const mutk::mutation::Model &model,
    size_t n, float t, Arg arg) {

    return model.CreateMatrix(n, t, arg);
}

template<typename Arg>
mutk::tensor_t create_transition_gamete_diploid_impl(const mutk::mutation::Model &model, 
    size_t n, float t, Arg arg) {

    using tensor_t = mutk::tensor_t;

    auto mat = create_transition_clone_haploid_impl(model, n, t, arg);
    tensor_t ret = tensor_t::from_shape({mutk::dim_width<1>(n),mutk::dim_width<2>(n)});

    // x/y -> a
    for(size_t u = 0; u < ret.size(); ++u) {
        auto s = xt::unravel_index(u, ret.shape());
        auto a = s[0];
        auto x = ALLELE[s[1]][0];
        auto y = ALLELE[s[1]][1];
        ret(u) = 0.5*(mat(a,x)+mat(a,y));
    }    
    return ret;  
}

template<int N, typename Arg>
inline
mutk::tensor_t create_transition_haploid_impl(const mutk::mutation::Model &model, 
    size_t n, float t, Arg arg) {
    if constexpr(N == 1) {
        return create_transition_clone_haploid_impl(model,n,t,arg);
    } else if constexpr(N == 2) {
        return create_transition_gamete_diploid_impl(model,n,t,arg);
    } else {
        static_assert(N == 1 || N == 2);
    }
}

struct clone_diploid_tag {};

template<int N, int M>
struct child_diploid_tag {};

template<int N>
struct child_selfing_tag {};

template<typename A>
struct transition_traits;

// G1/G2 => a/b
template<>
struct transition_traits<clone_diploid_tag> {
    static constexpr int MAT1_SIZE = 1;
    static constexpr int MAT2_SIZE = 1;
    static constexpr int TENSOR_SIZE = 2;

    static mutk::shape_t Dims(size_t n) {
        return {mutk::dim_width<2>(n), mutk::dim_width<2>(n)};
    }
    static size_t G1(const mutk::strides_t &coords) {
        return ALLELE[coords[1]][0];
    }
    static size_t G2(const mutk::strides_t &coords) {
        return ALLELE[coords[1]][1];
    }
};

// G1 x G2 => a/b
template<int N, int M>
struct transition_traits<child_diploid_tag<N, M>> {
    static constexpr int MAT1_SIZE = N;
    static constexpr int MAT2_SIZE = M;
    static constexpr int TENSOR_SIZE = 3;

    static mutk::shape_t  Dims(size_t n) {
        return {mutk::dim_width<2>(n), mutk::dim_width<N>(n), mutk::dim_width<M>(n)};
    }
    static size_t G1(const mutk::strides_t &coords) {
        return coords[1];
    }
    static size_t G2(const mutk::strides_t &coords) {
        return coords[2];
    }
};

// G1 x G1 => a/b
template<int N>
struct transition_traits<child_selfing_tag<N>> {
    static constexpr int MAT1_SIZE = N;
    static constexpr int MAT2_SIZE = N;
    static constexpr int TENSOR_SIZE = 2;

    static mutk::shape_t Dims(size_t n) {
        return {mutk::dim_width<2>(n), mutk::dim_width<N>(n)};
    }
    static size_t G1(const mutk::strides_t &coords) {
        return coords[1];
    }
    static size_t G2(const mutk::strides_t &coords) {
        return coords[1];
    }
};

template<typename Tag>
struct generate_transition_diploid {
    using traits = transition_traits<Tag>;

    template<typename Arg1, typename Arg2>
    inline
    mutk::tensor_t
    operator()(const mutk::mutation::Model &model, 
        size_t n, float t1, float t2, Arg1 arg1, Arg2 arg2) const {
     
        auto mat1 = create_transition_haploid_impl<traits::MAT1_SIZE>(model, n, t1, arg1);
        auto mat2 = create_transition_haploid_impl<traits::MAT2_SIZE>(model, n, t2, arg2);

        // G1/G2 => a/b OR G1 x G2 => a/b
        mutk::tensor_t ret = mutk::tensor_t::from_shape(traits::Dims(n));
        for(size_t u = 0; u < ret.size(); ++u) {
            auto s = xt::unravel_index(u, ret.shape());
            size_t a = ALLELE[s[0]][0];
            size_t b = ALLELE[s[0]][1];
            size_t g1 = traits::G1(s);
            size_t g2 = traits::G2(s);

            float f = mat1(a,g1) * mat2(b,g2);
            if(a != b) {
                f += mat1(b,g1) * mat2(a,g2);
            }
            ret(u) = f;
        }
        return ret;
    }
};

template<typename Tag, typename Arg>
auto create_transition_diploid_impl(const mutk::mutation::Model &model, 
    size_t n, float t1, float t2, Arg arg) {
    generate_transition_diploid<Tag> gen;

    if constexpr(std::is_same_v<Arg, mutk::mutation::any_t>) {
        return gen(model, n, t1, t2, arg, arg);
    } else if constexpr(std::is_same_v<Arg, mutk::mutation::mean_t>) {
        auto ret = gen(model, n, t1, t2, arg, mutk::mutation::ANY);
        ret += gen(model, n, t1, t2, mutk::mutation::ANY, arg);
        return ret;
    } else if constexpr(std::is_integral_v<Arg>) {
        auto ret = gen(model, n, t1, t2, 0, arg);
        for(Arg a=1; a <= arg; ++a) {
            ret += gen(model, n, t1, t2, a, arg-a);
        }
        return ret;
    } else {
        static_assert(std::is_integral_v<Arg>
            || std::is_same_v<Arg, mutk::mutation::mean_t>
            || std::is_same_v<Arg, mutk::mutation::any_t> );
    }     
}

template<int N, int M, typename Arg>
inline
mutk::tensor_t create_transition_child(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 2);
    assert(potential.axes.size() == 3);
    auto ret = create_transition_diploid_impl<child_diploid_tag<N,M>>(model, n, potential.parents[0].second,
        potential.parents[1].second, arg);
    return xt::transpose(ret, potential.axes);
}

template<int N, typename Arg>
inline
mutk::tensor_t create_transition_child_selfing(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    assert(potential.axes.size() == 2);
    auto ret = create_transition_diploid_impl<child_selfing_tag<N>>(model, n, potential.parents[0].second,
        potential.parents[0].second, arg);
    return xt::transpose(ret, potential.axes);
}

template<typename Arg>
mutk::tensor_t create_transition_clone_diploid(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    assert(potential.axes.size() == 2);
    auto ret = create_transition_diploid_impl<clone_diploid_tag>(model, n, potential.parents[0].second,
        potential.parents[0].second, arg);
    return xt::transpose(ret, potential.axes);
}

template<typename Arg>
mutk::tensor_t create_transition_clone_haploid(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    assert(potential.axes.size() == 2);
    auto ret =  create_transition_haploid_impl<1>(model, n, potential.parents[0].second, arg);
    return xt::transpose(ret, potential.axes);
}

template<typename Arg>
mutk::tensor_t create_transition_gamete_diploid(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    auto ret = create_transition_haploid_impl<2>(model, n, potential.parents[0].second, arg);
    return xt::transpose(ret, potential.axes);
}

template<typename Arg>
mutk::tensor_t create_transition_child_diploid_diploid(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 2);
    // don't shuffle here because the subfunction does it for us.
    return create_transition_child<2, 2>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_diploid_diploid") {
    auto test = [&](size_t n, float u, float k, std::vector<int> axes) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::PotentialType:: ChildDiploidDiploid, IndyId(0), IndyId(1), u, IndyId(2), 1.1f*u};
        pot.axes = axes;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_diploid_diploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);
        auto mat2 = model.CreateMatrix(n, 1.1f*u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 3);
        mutk::shape_t pot_shape = {n*(n+1)/2, n*(n+1)/2, n*(n+1)/2};
        REQUIRE(obs.shape(0) == pot_shape[pot.axes[0]]);
        REQUIRE(obs.shape(1) == pot_shape[pot.axes[1]]);
        REQUIRE(obs.shape(2) == pot_shape[pot.axes[2]]);
        for(int a1=0,a=0;a1<n;++a1) {
            for(int a2=0;a2<=a1;++a2,++a) {
                for(int b1=0,b=0;b1<n;++b1) {
                    for(int b2=0;b2<=b1;++b2,++b) {
                        for(int c1=0,c=0;c1<n;++c1) {
                            for(int c2=0;c2<=c1;++c2,++c) {
                                CAPTURE(a1);
                                CAPTURE(a2);
                                CAPTURE(a);
                                CAPTURE(b1);
                                CAPTURE(b2);
                                CAPTURE(b);
                                CAPTURE(c1);
                                CAPTURE(c2);
                                CAPTURE(c);
                                // b1/b2 x c1/c2 -> a1/a2
                                float expected1 = 0.5f*mat1(a1,b1) + 0.5f*mat1(a1,b2);
                                float expected2 = 0.5f*mat2(a2,c1) + 0.5f*mat2(a2,c2);
                                float expected = expected1*expected2;
                                if(a1 != a2) {
                                    expected1 = 0.5f*mat1(a2,b1) + 0.5f*mat1(a2,b2);
                                    expected2 = 0.5f*mat2(a1,c1) + 0.5f*mat2(a1,c2);
                                    expected += expected1*expected2;
                                }

                                std::array<int, 3> indexes = {a,b,c};
                                int o1 = indexes[pot.axes[0]];
                                int o2 = indexes[pot.axes[1]];
                                int o3 = indexes[pot.axes[2]];
                                CHECK(obs(o1,o2,o3) == doctest::Approx(expected));
                            }
                        }
                    }
                }
            }
        }
    };
    SUBCASE("Axes are {0, 1, 2}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {0, 1, 2});
        };
        run_mutation_tests(subtest);
    }
    SUBCASE("Axes are {1, 0, 2}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {1, 0, 2});
        };
        run_mutation_tests(subtest);
    }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_child_haploid_diploid(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 2);
    // no transpose here
    return create_transition_child<1, 2>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_haploid_diploid") {
    auto test = [&](size_t n, float u, float k, std::vector<int> axes) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::PotentialType:: ChildHaploidDiploid, IndyId(0), IndyId(1), u, IndyId(2), 1.1f*u};
        pot.axes = axes;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_haploid_diploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);
        auto mat2 = model.CreateMatrix(n, 1.1f*u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 3);
        mutk::shape_t pot_shape = {n*(n+1)/2, n, n*(n+1)/2};
        REQUIRE(obs.shape(0) == pot_shape[pot.axes[0]]);
        REQUIRE(obs.shape(1) == pot_shape[pot.axes[1]]);
        REQUIRE(obs.shape(2) == pot_shape[pot.axes[2]]);
        for(int a1=0,a=0;a1<n;++a1) {
            for(int a2=0;a2<=a1;++a2,++a) {
                for(int b=0;b<n;++b) {
                    for(int c1=0,c=0;c1<n;++c1) {
                        for(int c2=0;c2<=c1;++c2,++c) {
                            CAPTURE(a1);
                            CAPTURE(a2);
                            CAPTURE(a);
                            CAPTURE(b);
                            CAPTURE(c1);
                            CAPTURE(c2);
                            CAPTURE(c);
                            // b x c1/c2 -> a1/a2
                            float expected1 = mat1(a1,b);
                            float expected2 = 0.5f*mat2(a2,c1) + 0.5f*mat2(a2,c2);
                            float expected = expected1*expected2;
                            if(a1 != a2) {
                                expected1 = mat1(a2,b);
                                expected2 = 0.5f*mat2(a1,c1) + 0.5f*mat2(a1,c2);
                                expected += expected1*expected2;
                            }

                            std::array<int, 3> indexes = {a,b,c};
                            int o1 = indexes[pot.axes[0]];
                            int o2 = indexes[pot.axes[1]];
                            int o3 = indexes[pot.axes[2]];
                            CHECK(obs(o1,o2,o3) == doctest::Approx(expected));
                        }
                    }
                }
            }
        }
    };
    SUBCASE("Axes are {0, 1, 2}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {0, 1, 2});
        };
        run_mutation_tests(subtest);
    }
    SUBCASE("Axes are {1, 0, 2}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {1, 0, 2});
        };
        run_mutation_tests(subtest);
    }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_child_diploid_haploid(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 2);
    // no transpose here
    return create_transition_child<2, 1>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_diploid_haploid") {
    auto test = [&](size_t n, float u, float k, std::vector<int> axes) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::PotentialType:: ChildDiploidHaploid, IndyId(0), IndyId(1), u, IndyId(2), 1.1f*u};
        pot.axes = axes;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_diploid_haploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);
        auto mat2 = model.CreateMatrix(n, 1.1f*u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 3);
        mutk::shape_t pot_shape = {n*(n+1)/2, n*(n+1)/2, n};
        REQUIRE(obs.shape(0) == pot_shape[pot.axes[0]]);
        REQUIRE(obs.shape(1) == pot_shape[pot.axes[1]]);
        REQUIRE(obs.shape(2) == pot_shape[pot.axes[2]]);
        for(int a1=0,a=0;a1<n;++a1) {
            for(int a2=0;a2<=a1;++a2,++a) {
                for(int b1=0,b=0;b1<n;++b1) {
                    for(int b2=0;b2<=b1;++b2,++b) {
                        for(int c=0;c<n;++c) {
                            CAPTURE(a1);
                            CAPTURE(a2);
                            CAPTURE(a);
                            CAPTURE(b1);
                            CAPTURE(b2);
                            CAPTURE(b);
                            CAPTURE(c);
                            // b1/b2 x c -> a1/a2
                            float expected1 = 0.5f*mat1(a1,b1) + 0.5f*mat1(a1,b2);
                            float expected2 = mat2(a2,c);
                            float expected = expected1*expected2;
                            if(a1 != a2) {
                                expected1 = 0.5f*mat1(a2,b1) + 0.5f*mat1(a2,b2);
                                expected2 = mat2(a1,c);
                                expected += expected1*expected2;
                            }

                            std::array<int, 3> indexes = {a,b,c};
                            int o1 = indexes[pot.axes[0]];
                            int o2 = indexes[pot.axes[1]];
                            int o3 = indexes[pot.axes[2]];
                            CHECK(obs(o1,o2,o3) == doctest::Approx(expected));
                        }
                    }
                }
            }
        }
    };
    SUBCASE("Axes are {0, 1, 2}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {0, 1, 2});
        };
        run_mutation_tests(subtest);
    }
    SUBCASE("Axes are {1, 0, 2}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {1, 0, 2});
        };
        run_mutation_tests(subtest);
    }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_child_haploid_haploid(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 2);
    // no transpose here
    return create_transition_child<1, 1>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_haploid_haploid") {
    auto test = [&](size_t n, float u, float k, std::vector<int> axes) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::PotentialType:: ChildHaploidHaploid, IndyId(0), IndyId(1), u, IndyId(2), 1.1f*u};
        pot.axes = axes;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_haploid_haploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);
        auto mat2 = model.CreateMatrix(n, 1.1f*u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 3);
        mutk::shape_t pot_shape = {n*(n+1)/2, n, n};
        REQUIRE(obs.shape(0) == pot_shape[pot.axes[0]]);
        REQUIRE(obs.shape(1) == pot_shape[pot.axes[1]]);
        REQUIRE(obs.shape(2) == pot_shape[pot.axes[2]]);
        for(int a1=0,a=0;a1<n;++a1) {
            for(int a2=0;a2<=a1;++a2,++a) {
                for(int b=0;b<n;++b) {
                    for(int c=0;c<n;++c) {
                        CAPTURE(a1);
                        CAPTURE(a2);
                        CAPTURE(a);
                        CAPTURE(b);
                        CAPTURE(c);
                        // b x c -> a1/a2
                        float expected1 = mat1(a1,b);
                        float expected2 = mat2(a2,c);
                        float expected = expected1*expected2;
                        if(a1 != a2) {
                            expected1 = mat1(a2,b);
                            expected2 = mat2(a1,c);
                            expected += expected1*expected2;
                        }

                        std::array<int, 3> indexes = {a,b,c};
                        int o1 = indexes[pot.axes[0]];
                        int o2 = indexes[pot.axes[1]];
                        int o3 = indexes[pot.axes[2]];
                        CHECK(obs(o1,o2,o3) == doctest::Approx(expected));
                    }
                }
            }
        }
    };
    SUBCASE("Axes are {0, 1, 2}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {0, 1, 2});
        };
        run_mutation_tests(subtest);
    }
    SUBCASE("Axes are {1, 0, 2}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {1, 0, 2});
        };
        run_mutation_tests(subtest);
    }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_child_selfing_diploid(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    // no transpose here
    return create_transition_child_selfing<2>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_selfing_diploid") {
    auto test = [&](size_t n, float u, float k, std::vector<int> axes) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::PotentialType::ChildSelfingDiploid, IndyId(0), IndyId(1), u};
        pot.axes = axes;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_selfing_diploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 2);
        mutk::shape_t pot_shape = {n*(n+1)/2, n*(n+1)/2};
        REQUIRE(obs.shape(0) == pot_shape[pot.axes[0]]);
        REQUIRE(obs.shape(1) == pot_shape[pot.axes[1]]);
        for(int a1=0,a=0;a1<n;++a1) {
            for(int a2=0;a2<=a1;++a2,++a) {
                for(int b1=0,b=0;b1<n;++b1) {
                    for(int b2=0;b2<=b1;++b2,++b) {
                        CAPTURE(a1);
                        CAPTURE(a2);
                        CAPTURE(a);
                        CAPTURE(b1);
                        CAPTURE(b2);
                        CAPTURE(b);
                        // b1/b2 x b1/b2 -> a1/a2
                        float expected1 = 0.5f*mat1(a1,b1) + 0.5f*mat1(a1,b2);
                        float expected2 = 0.5f*mat1(a2,b1) + 0.5f*mat1(a2,b2);
                        float expected = expected1*expected2;
                        if(a1 != a2) {
                            expected1 = 0.5f*mat1(a2,b1) + 0.5f*mat1(a2,b2);
                            expected2 = 0.5f*mat1(a1,b1) + 0.5f*mat1(a1,b2);
                            expected += expected1*expected2;
                        }

                        std::array<int, 2> indexes = {a,b};
                        int o1 = indexes[pot.axes[0]];
                        int o2 = indexes[pot.axes[1]];
                        CHECK(obs(o1,o2) == doctest::Approx(expected));
                    }
                }
            }
        }
    };
    SUBCASE("Axes are {0, 1}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {0, 1});
        };
        run_mutation_tests(subtest);
    }
    SUBCASE("Axes are {1, 0}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {1, 0});
        };
        run_mutation_tests(subtest);
    }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_child_selfing_haploid(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    // no transpose here
    return create_transition_child_selfing<1>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_selfing_haploid") {
    auto test = [&](size_t n, float u, float k, std::vector<int> axes) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::PotentialType::ChildSelfingHaploid, IndyId(0), IndyId(1), u};
        pot.axes = axes;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_selfing_haploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 2);
        mutk::shape_t pot_shape = {n*(n+1)/2, n};
        REQUIRE(obs.shape(0) == pot_shape[pot.axes[0]]);
        REQUIRE(obs.shape(1) == pot_shape[pot.axes[1]]);
        for(int a1=0,a=0;a1<n;++a1) {
            for(int a2=0;a2<=a1;++a2,++a) {
                for(int b=0;b<n;++b) {
                    CAPTURE(a1);
                    CAPTURE(a2);
                    CAPTURE(a);
                    CAPTURE(b);
                    // b x b -> a1/a2
                    float expected1 = mat1(a1,b);
                    float expected2 = mat1(a2,b);
                    float expected = expected1*expected2;
                    if(a1 != a2) {
                        expected1 = mat1(a2,b);
                        expected2 = mat1(a1,b);
                        expected += expected1*expected2;
                    }

                    std::array<int, 2> indexes = {a,b};
                    int o1 = indexes[pot.axes[0]];
                    int o2 = indexes[pot.axes[1]];
                    CHECK(obs(o1,o2) == doctest::Approx(expected));
                }
            }
        }
    };
    SUBCASE("Axes are {0, 1}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {0, 1});
        };
        run_mutation_tests(subtest);
    }
    SUBCASE("Axes are {1, 0}.") {
        auto subtest = [&](size_t n, float u, float k) {
            return test(n, u, k, {1, 0});
        };
        run_mutation_tests(subtest);
    }
}
// LCOV_EXCL_STOP

namespace mutk {
template<typename Arg>
inline
mutk::tensor_t create_mutation_potential(const mutk::mutation::Model &model, 
    size_t n, const potential_t &potential, Arg arg) {
    using P = mutk::PotentialType;

    switch(potential.type) {
    case P::CloneDiploid:
        return create_transition_clone_diploid(model, n, potential, arg);
    case P::CloneHaploid:
        return create_transition_clone_haploid(model, n, potential, arg);
    case P::GameteDiploid:
        return create_transition_gamete_diploid(model, n, potential, arg);
    case P::ChildDiploidDiploid:
        return create_transition_child_diploid_diploid(model, n, potential, arg);
    case P::ChildHaploidDiploid:
        return create_transition_child_haploid_diploid(model, n, potential, arg);
    case P::ChildDiploidHaploid:
        return create_transition_child_diploid_haploid(model, n, potential, arg);
    case P::ChildHaploidHaploid:
        return create_transition_child_haploid_haploid(model, n, potential, arg);
    case P::ChildSelfingDiploid:
        return create_transition_child_selfing_diploid(model, n, potential, arg);
    case P::ChildSelfingHaploid:
        return create_transition_child_selfing_haploid(model, n, potential, arg);
    case P::Unit:
        return create_unit_potential(n, potential);
    case P::FounderDiploid:
        return model.CreatePriorDiploid(n);
    case P::FounderHaploid: 
        return model.CreatePriorHaploid(n);
    default:
        break;
    };
    return {};
}
}

mutk::tensor_t mutk::mutation::Model::CreatePotential(size_t n, const potential_t &potential, any_t arg) {
    return create_mutation_potential(*this, n, potential, arg);
}

mutk::tensor_t mutk::mutation::Model::CreatePotential(size_t n, const potential_t &potential, mean_t arg) {
    return create_mutation_potential(*this, n, potential, arg);
}

mutk::tensor_t mutk::mutation::Model::CreatePotential(size_t n, const potential_t &potential, size_t arg) {
    return create_mutation_potential(*this, n, potential, arg);
}

#endif
