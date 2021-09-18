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
#include <doctest/doctest.h>

#include <mutk/mutation.hpp>
#include <mutk/relationship_graph.hpp>

// Libraries needed for testing
#include <boost/numeric/ublas/matrix.hpp>

namespace {
// Structure useful for unit testing
struct kalleles_test_mat {
    using mat_t = boost::numeric::ublas::matrix<float,
        boost::numeric::ublas::column_major,std::vector<float>>;

    mat_t Q;

    kalleles_test_mat(int n, float k) : Q(n+1,n+1) {
        using namespace boost::numeric::ublas;
        // Build matrix
        vector<float> freqs(n+1, 1.0/k);
        freqs[n] = 1.0-n/k;

        for(int i=0; i<=n; ++i) {
            for(int j=0; j<=n; ++j) {
                Q(i,j) = freqs[j];
            }
            Q(i,i) = Q(i,i)-1.0;
        }

        // Scale matrix into substitution time
        Q *= k/(k-1);
    }
};

template<typename F>
void run_mutation_tests(F test) {
    test(2, 0.0,  4.0);
    test(4, 0.0,  4.0);

    test(1, 1e-8, 4.0);
    test(2, 1e-8, 4.0);
    test(3, 1e-8, 4.0);
    test(4, 1e-8, 4.0);
    
    test(4, 1e-8, 4.0);
    test(4, 1e-9, 5.0);
    test(4, 1e-6, 6.0);
    test(4, 1e-3, 7.0);
}

template<typename V>
inline int shuffle_pos(int i, const V & v) {
    auto it = std::find(std::begin(v), std::end(v), i);
    return std::distance(std::begin(v), it);
}

} // anon namespace

namespace {
constexpr int ALLELE[][2] = {
    {0,0},
    {0,1}, {1,1},
    {0,2}, {1,2}, {2,2},
    {0,3}, {1,3}, {2,3}, {3,3},
    {0,4}, {1,4}, {2,4}, {3,4}, {4,4},
};
}

using KAllelesModel = mutk::mutation::KAllelesModel;
using potential_t = mutk::RelationshipGraph::potential_t;

// LCOV_EXCL_START
TEST_CASE("KAllelesModel-Constructor") {
    CHECK_NOTHROW(KAllelesModel(4.0, 0.001, 0.0, 0.0, 0.0));
    CHECK_THROWS_AS(KAllelesModel(1.0, 0.001, 0.0, 0.0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(KAllelesModel(4.0, -0.1, 0.0, 0.0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(KAllelesModel(4.0, 0.001, 1.1, 0.0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(KAllelesModel(4.0, 0.001, -3000, 0.0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(KAllelesModel(4.0, 0.001, 0.0, 1.1, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(KAllelesModel(4.0, 0.001, 0.0, -4000, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(KAllelesModel(4.0, 0.001, 0.0, 0.0, 1.1), std::invalid_argument);
    CHECK_THROWS_AS(KAllelesModel(4.0, 0.001, 0.0, 0.0, -3000), std::invalid_argument);
}
// LCOV_EXCL_STOP

// ret(i,j) = P(i|j)
KAllelesModel::tensor_t KAllelesModel::CreateMatrix(size_t n, float t, any_t) const {
    assert(n > 0);
    assert(n <= 5);

    double beta = t*k_/(k_-1.0);
    double p_ji = -1.0/k_*expm1(-beta);
    double p_jj = exp(-beta) + p_ji;

    tensor_t ret = tensor_t::from_shape({n,n});

    for(size_t i = 0; i < n; ++i) {
        for(size_t j = 0; j < n; ++j) {
            ret(i,j) = (i == j) ? p_jj : p_ji;
        }
    }
    return ret;
}

// LCOV_EXCL_START
TEST_CASE("KAllelesModel-CreateMatrix with any_t") {
    using namespace boost::numeric::ublas;
    
    auto test = [&](int n, float u, float k) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        KAllelesModel model(k, 0.001, 0, 0, 0);

        using mat_t = kalleles_test_mat::mat_t;

        kalleles_test_mat mat(n,k);

        mat_t P = identity_matrix<float>(n+1);
        mat_t m = P;
        float t = 1.0;
        float f = 1.0;
        // Use the first 11 elements of the expm series 
        for(int g=1;g <= 11; ++g) {
            mat_t a = prec_prod(m,mat.Q);
            m = a;
            t *= u;
            f *= g;
            P += m * (t/f);
        }
        auto obs = model.CreateMatrix(n, u, mutk::mutation::ANY);
        REQUIRE(obs.dimension() == 2);
        REQUIRE(obs.shape(0) == n);
        REQUIRE(obs.shape(1) == n);
        for(int i=0; i < n; ++i) {
            for(int j=0; j < n; ++j) {
                CAPTURE(i);
                CAPTURE(j);
                CHECK(obs(i,j) == doctest::Approx(P(i,j)));
            }
        }
    };
    run_mutation_tests(test);
}
// LCOV_EXCL_STOP

// ret(i,j) = P(i & x mutations | j)
KAllelesModel::tensor_t KAllelesModel::CreateMatrix(size_t n, float t, size_t x) const {
    assert(n > 0);
    assert(n <= 5);
    assert(x >= 0);
    
    double p_x;
    if(t == 0.0) {
        p_x = (x==0) ? 1.0 : 0.0;
    } else {
        p_x = exp(-t+x*log(t)-lgamma(x+1));
    }

    double p_ji = (1.0-pow(-1.0/(k_-1.0),x))/k_;
    double p_jj = (1.0+(k_-1.0)*pow(-1.0/(k_-1.0),x))/k_;

    tensor_t ret = tensor_t::from_shape({n,n});

    for(size_t i = 0; i < n; ++i) {
        for(size_t j = 0; j < n; ++j) {
            ret(i,j) = (i == j) ? p_x*p_jj : p_x*p_ji;
        }
    }
    return ret;
}

// LCOV_EXCL_START
TEST_CASE("KAllelesModel-CreateMatrix with int") {
    using namespace boost::numeric::ublas;
    
    auto test = [&](int n, float u, float k) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        KAllelesModel model(k, 0.001, 0, 0, 0);

        using mat_t = kalleles_test_mat::mat_t;

        kalleles_test_mat mat(n,k);

        mat_t P = identity_matrix<float>(n+1);
        mat_t m = P;
        mat_t J = mat.Q+m;
        float t;
        float f;
        for(int x=0; x <= 10; ++x) {
            CAPTURE(x);
            if(x == 0) {
                P = m*exp(-u);
                f = 1.0;
                t = 1.0;
            } else {
                mat_t a = prec_prod(m,J);
                m = a;
                t *= u;
                f *= x;
                P = m*(exp(-u)*t/f);
            }
            auto obs = model.CreateMatrix(n, u, x);
            REQUIRE(obs.dimension() == 2);
            REQUIRE(obs.shape(0) == n);
            REQUIRE(obs.shape(1) == n);
            for(int i=0; i < n; ++i) {
                for(int j=0; j < n; ++j) {
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

// ret(i,j) = E[num of mutations | i,j]*P(i|j)
KAllelesModel::tensor_t KAllelesModel::CreateMatrix(size_t n, float t, mean_t) const {
    assert(n > 0);
    assert(n <= 5);

    double beta = k_*t/(k_-1.0);
    double p_jj = -t/k_*expm1(-beta);
    double p_ji = (t-p_jj)/(k_-1.0);

    tensor_t ret = tensor_t::from_shape({n,n});

    for(size_t i = 0; i < n; ++i) {
        for(size_t j = 0; j < n; ++j) {
            ret(i,j) = (i == j) ? p_jj : p_ji;
        }
    }
    return ret;
}

// LCOV_EXCL_START
TEST_CASE("KAllelesModel-CreateMatrix with mean_t") {
    using namespace boost::numeric::ublas;
    
    auto test = [&](int n, float u, float k) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        KAllelesModel model(k, 0.001, 0, 0, 0);

        using mat_t = kalleles_test_mat::mat_t;

        kalleles_test_mat mat(n,k);

        mat_t P = identity_matrix<float>(n+1);
        mat_t m = P;
        mat_t J = mat.Q+m;
        mat_t S;
        float t;
        float f;
        for(int x=0; x <= 10; ++x) {
            if(x == 0) {
                P = m*exp(-u);
                S = P*x;
                f = 1.0;
                t = 1.0;
            } else {
                mat_t a = prec_prod(m,J);
                m = a;
                t *= u;
                f *= x;
                P = m*(exp(-u)*t/f);
                S += P*x;
            }
        }
        auto obs = model.CreateMatrix(n, u, mutk::mutation::MEAN);
        REQUIRE(obs.dimension() == 2);
        REQUIRE(obs.shape(0) == n);
        REQUIRE(obs.shape(1) == n);
        for(int i=0; i < n; ++i) {
            for(int j=0; j < n; ++j) {
                CAPTURE(i);
                CAPTURE(j);
                CHECK(obs(i,j) == doctest::Approx(S(i,j)));
            }
        }
    };

    run_mutation_tests(test);
}
// LCOV_EXCL_STOP

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
    auto test_haploid = [](int n, float theta, float hap_bias,
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
    auto test_diploid = [](int n, float theta, float hom_bias,
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
    int n, float t, Arg arg) {

    return model.CreateMatrix(n, t, arg);
}

template<typename Arg>
mutk::tensor_t create_transition_gamete_diploid_impl(const mutk::mutation::Model &model, 
    int n, float t, Arg arg) {

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
    int n, float t1, float t2, Arg arg) {
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
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 2);
    assert(potential.shuffle.size() == 3);
    //auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1],potential.shuffle[2]);
    return create_transition_diploid_impl<child_diploid_tag<N,M>>(model, n, potential.parents[0].second,
        potential.parents[1].second, arg);
}

template<int N, typename Arg>
inline
mutk::tensor_t create_transition_child_selfing(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    assert(potential.shuffle.size() == 2);
    //auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);
    return create_transition_diploid_impl<child_selfing_tag<N>>(model, n, potential.parents[0].second,
        potential.parents[0].second, arg);
}

template<typename Arg>
mutk::tensor_t create_transition_clone_diploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    assert(potential.shuffle.size() == 2);
    //auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);
    return create_transition_diploid_impl<clone_diploid_tag>(model, n, potential.parents[0].second,
        potential.parents[0].second, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_clone_diploid") {
    SUBCASE("Passing mutk::mutation::ANY as argument.") {
        auto test = [&](int n, float u, float k, std::vector<int> shuf) {
            CAPTURE(n);
            CAPTURE(k);
            CAPTURE(u);

            potential_t pot{mutk::detail::Potential::CloneDiploid, 0, 1, u};
            pot.shuffle = shuf;
            KAllelesModel model{k, 0.001, 0, 0, 0};
            auto obs = create_transition_clone_diploid(model, n, pot, mutk::mutation::ANY);

            auto mat = model.CreateMatrix(n, u, mutk::mutation::ANY);

            REQUIRE(obs.dimension() == 2);
            REQUIRE(obs.shape(shuffle_pos(0,pot.shuffle)) == n*(n+1)/2);
            REQUIRE(obs.shape(shuffle_pos(1,pot.shuffle)) == n*(n+1)/2);
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
                            // b1/b2 -> a1/a2
                            float expected = mat(a1,b1)*mat(a2,b2);
                            if(a1 != a2) {
                                expected += mat(a2,b1)*mat(a1,b2);
                            }
                            std::array<int, 2> indexes = {a,b};
                            int o1 = indexes[pot.shuffle[0]];
                            int o2 = indexes[pot.shuffle[1]];
                            CHECK(obs(o1,o2) == doctest::Approx(expected));
                        }
                    }
                }
            }
        };
        SUBCASE("Shuffle order is {0, 1}.") {
            auto subtest = [&](int n, float u, float k) {
                return test(n, u, k, {0,1});
            };
            run_mutation_tests(subtest);
        }
        // SUBCASE("Shuffle order is {1, 0}.") {
        //     auto subtest = [&](int n, float u, float k) {
        //         return test(n, u, k, {1,0});
        //     };
        //     run_mutation_tests(subtest);
        // }
    }
    SUBCASE("Passing mutk::mutation::MEAN as argument.") {
        auto test = [&](int n, float u, float k, std::vector<int> shuf) {
            CAPTURE(n);
            CAPTURE(k);
            CAPTURE(u);

            potential_t pot{mutk::detail::Potential::CloneDiploid, 0, 1, u};
            pot.shuffle = shuf;
            KAllelesModel model{k, 0.001, 0, 0, 0};
            auto obs = create_transition_clone_diploid(model, n, pot, mutk::mutation::MEAN);

            auto mat1 = model.CreateMatrix(n, u, mutk::mutation::MEAN);
            auto mat2 = model.CreateMatrix(n, u, mutk::mutation::ANY);

            REQUIRE(obs.dimension() == 2);
            REQUIRE(obs.shape(shuffle_pos(0,pot.shuffle)) == n*(n+1)/2);
            REQUIRE(obs.shape(shuffle_pos(1,pot.shuffle)) == n*(n+1)/2);
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
                            // b1/b2 -> a1/a2
                            float expected = mat1(a1,b1)*mat2(a2,b2) + mat2(a1,b1)*mat1(a2,b2);
                            if(a1 != a2) {
                                expected += mat1(a2,b1)*mat2(a1,b2) + mat2(a2,b1)*mat1(a1,b2);
                            }
                            std::array<int, 2> indexes = {a,b};
                            int o1 = indexes[pot.shuffle[0]];
                            int o2 = indexes[pot.shuffle[1]];
                            CHECK(obs(o1,o2) == doctest::Approx(expected));
                        }
                    }
                }
            }
        };
        SUBCASE("Shuffle order is {0, 1}.") {
            auto subtest = [&](int n, float u, float k) {
                return test(n, u, k, {0,1});
            };
            run_mutation_tests(subtest);
        }
        // SUBCASE("Shuffle order is {1, 0}.") {
        //     auto subtest = [&](int n, float u, float k) {
        //         return test(n, u, k, {1,0});
        //     };
        //     run_mutation_tests(subtest);
        // }
    }
    SUBCASE("Passing an integer as argument.") {
        auto test = [&](int n, float u, float k, std::vector<int> shuf, int val) {
            CAPTURE(n);
            CAPTURE(k);
            CAPTURE(u);

            potential_t pot{mutk::detail::Potential::CloneDiploid, 0, 1, u};
            pot.shuffle = shuf;
            KAllelesModel model{k, 0.001, 0, 0, 0};
            auto obs = create_transition_clone_diploid(model, n, pot, val);

            REQUIRE(obs.dimension() == 2);
            REQUIRE(obs.shape(shuffle_pos(0,pot.shuffle)) == n*(n+1)/2);
            REQUIRE(obs.shape(shuffle_pos(1,pot.shuffle)) == n*(n+1)/2);
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
                            // b1/b2 -> a1/a2
                            float expected = 0.0;
                            for(int i=0;i<=val;++i) {           
                                auto mat1 = model.CreateMatrix(n, u, i);
                                auto mat2 = model.CreateMatrix(n, u, val-i);
                                expected += mat1(a1,b1)*mat2(a2,b2);
                                if(a1 != a2) {
                                    expected += mat1(a2,b1)*mat2(a1,b2);
                                }
                            }
                            std::array<int, 2> indexes = {a,b};
                            int o1 = indexes[pot.shuffle[0]];
                            int o2 = indexes[pot.shuffle[1]];
                            CHECK(obs(o1,o2) == doctest::Approx(expected));
                        }
                    }
                }
            }
        };
        SUBCASE("Shuffle order is {0, 1}. Integer is 0.") {
            auto subtest = [&](int n, float u, float k) {
                return test(n, u, k, {0,1}, 0);
            };
            run_mutation_tests(subtest);
        }
        // SUBCASE("Shuffle order is {1, 0}. Integer is 1.") {
        //     auto subtest = [&](int n, float u, float k) {
        //         return test(n, u, k, {1,0}, 1);
        //     };
        //     run_mutation_tests(subtest);
        // }
    }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_clone_haploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    assert(potential.shuffle.size() == 2);
    //auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);
    return create_transition_haploid_impl<1>(model, n, potential.parents[0].second, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_clone_haploid") {
    auto test = [&](int n, float u, float k, std::vector<int> shuf) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::detail::Potential::CloneHaploid, 0, 1, u};
        pot.shuffle = shuf;
        KAllelesModel model{k, 0.001, 0, 0, 0};
        auto obs = create_transition_clone_haploid(model, n, pot, mutk::mutation::ANY);

        auto mat = model.CreateMatrix(n, u, mutk::mutation::ANY);
        REQUIRE(obs.dimension() == 2);
        REQUIRE(obs.shape(shuffle_pos(0,pot.shuffle)) == n);
        REQUIRE(obs.shape(shuffle_pos(1,pot.shuffle)) == n);

        for(int a=0;a<n;++a) {
            for(int b=0;b<n;++b) {
                CAPTURE(a);
                CAPTURE(b);
                // b -> a
                float expected = mat(a,b);
                std::array<int, 2> indexes = {a,b};
                int o1 = indexes[pot.shuffle[0]];
                int o2 = indexes[pot.shuffle[1]];
                CHECK(obs(o1,o2) == doctest::Approx(expected));
            }
        }
    };
    SUBCASE("Shuffle order is {0, 1}.") {
        auto subtest = [&](int n, float u, float k) {
            return test(n, u, k, {0,1});
        };
        run_mutation_tests(subtest);
    }
    // SUBCASE("Shuffle order is {1, 0}.") {
    //     auto subtest = [&](int n, float u, float k) {
    //         return test(n, u, k, {1,0});
    //     };
    //     run_mutation_tests(subtest);
    // }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_gamete_diploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    //auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);
    return create_transition_haploid_impl<2>(model, n, potential.parents[0].second, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_gamete_diploid") {
    auto test = [&](int n, float u, float k, std::vector<int> shuf) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::detail::Potential::GameteDiploid, 0, 1, u};
        pot.shuffle = shuf;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_gamete_diploid(model, n, pot, mutk::mutation::ANY);
        auto mat = model.CreateMatrix(n, u, mutk::mutation::ANY);
        REQUIRE(obs.dimension() == 2);
        REQUIRE(obs.shape(shuffle_pos(0,pot.shuffle)) == n);
        REQUIRE(obs.shape(shuffle_pos(1,pot.shuffle)) == n*(n+1)/2);        
        for(int a=0;a<n;++a) {
            for(int b1=0,b=0;b1<n;++b1) {
                for(int b2=0;b2<=b1;++b2,++b) {
                    CAPTURE(a);
                    CAPTURE(b);
                    CAPTURE(b1);
                    CAPTURE(b2);
                    // b1/b2 -> a
                    float expected = 0.5*mat(a,b1) + 0.5*mat(a,b2);
                    std::array<int, 2> indexes = {a,b};
                    int o1 = indexes[pot.shuffle[0]];
                    int o2 = indexes[pot.shuffle[1]];
                    CHECK(obs(o1,o2) == doctest::Approx(expected));
                }
            }
        }
    };
    SUBCASE("Shuffle order is {0, 1}.") {
        auto subtest = [&](int n, float u, float k) {
            return test(n, u, k, {0,1});
        };
        run_mutation_tests(subtest);
    }
    // SUBCASE("Shuffle order is {1, 0}.") {
    //     auto subtest = [&](int n, float u, float k) {
    //         return test(n, u, k, {1,0});
    //     };
    //     run_mutation_tests(subtest);
    // }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_child_diploid_diploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 2);
    //auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1],potential.shuffle[2]);
    return create_transition_child<2, 2>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_diploid_diploid") {
    auto test = [&](int n, float u, float k, std::vector<int> shuf) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::detail::Potential:: ChildDiploidDiploid, 0, 1, u, 2, 1.1f*u};
        pot.shuffle = shuf;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_diploid_diploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);
        auto mat2 = model.CreateMatrix(n, 1.1f*u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 3);
        REQUIRE(obs.shape(shuffle_pos(0,pot.shuffle)) == n*(n+1)/2);
        REQUIRE(obs.shape(shuffle_pos(1,pot.shuffle)) == n*(n+1)/2);
        REQUIRE(obs.shape(shuffle_pos(2,pot.shuffle)) == n*(n+1)/2);
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
                                int o1 = indexes[pot.shuffle[0]];
                                int o2 = indexes[pot.shuffle[1]];
                                int o3 = indexes[pot.shuffle[2]];
                                CHECK(obs(o1,o2,o3) == doctest::Approx(expected));
                            }
                        }
                    }
                }
            }
        }
    };
    SUBCASE("Shuffle order is {0, 1, 2}.") {
        auto subtest = [&](int n, float u, float k) {
            return test(n, u, k, {0, 1, 2});
        };
        run_mutation_tests(subtest);
    }
    // SUBCASE("Shuffle order is {1, 0, 2}.") {
    //     auto subtest = [&](int n, float u, float k) {
    //         return test(n, u, k, {1, 0, 2});
    //     };
    //     run_mutation_tests(subtest);
    // }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_child_haploid_diploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 2);
    //auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1],potential.shuffle[2]);
    return create_transition_child<1, 2>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_haploid_diploid") {
    auto test = [&](int n, float u, float k, std::vector<int> shuf) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::detail::Potential:: ChildHaploidDiploid, 0, 1, u, 2, 1.1f*u};
        pot.shuffle = shuf;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_haploid_diploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);
        auto mat2 = model.CreateMatrix(n, 1.1f*u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 3);
        REQUIRE(obs.shape(shuffle_pos(0,pot.shuffle)) == n*(n+1)/2);
        REQUIRE(obs.shape(shuffle_pos(1,pot.shuffle)) == n);
        REQUIRE(obs.shape(shuffle_pos(2,pot.shuffle)) == n*(n+1)/2);
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
                            int o1 = indexes[pot.shuffle[0]];
                            int o2 = indexes[pot.shuffle[1]];
                            int o3 = indexes[pot.shuffle[2]];
                            CHECK(obs(o1,o2,o3) == doctest::Approx(expected));
                        }
                    }
                }
            }
        }
    };
    SUBCASE("Shuffle order is {0, 1, 2}.") {
        auto subtest = [&](int n, float u, float k) {
            return test(n, u, k, {0, 1, 2});
        };
        run_mutation_tests(subtest);
    }
    // SUBCASE("Shuffle order is {1, 0, 2}.") {
    //     auto subtest = [&](int n, float u, float k) {
    //         return test(n, u, k, {1, 0, 2});
    //     };
    //     run_mutation_tests(subtest);
    // }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_child_diploid_haploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 2);
    //auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1],potential.shuffle[2]);
    return create_transition_child<2, 1>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_diploid_haploid") {
    auto test = [&](int n, float u, float k, std::vector<int> shuf) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::detail::Potential:: ChildDiploidHaploid, 0, 1, u, 2, 1.1f*u};
        pot.shuffle = shuf;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_diploid_haploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);
        auto mat2 = model.CreateMatrix(n, 1.1f*u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 3);
        REQUIRE(obs.shape(shuffle_pos(0,pot.shuffle)) == n*(n+1)/2);
        REQUIRE(obs.shape(shuffle_pos(1,pot.shuffle)) == n*(n+1)/2);
        REQUIRE(obs.shape(shuffle_pos(2,pot.shuffle)) == n);
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
                            int o1 = indexes[pot.shuffle[0]];
                            int o2 = indexes[pot.shuffle[1]];
                            int o3 = indexes[pot.shuffle[2]];
                            CHECK(obs(o1,o2,o3) == doctest::Approx(expected));
                        }
                    }
                }
            }
        }
    };
    SUBCASE("Shuffle order is {0, 1, 2}.") {
        auto subtest = [&](int n, float u, float k) {
            return test(n, u, k, {0, 1, 2});
        };
        run_mutation_tests(subtest);
    }
    // SUBCASE("Shuffle order is {1, 0, 2}.") {
    //     auto subtest = [&](int n, float u, float k) {
    //         return test(n, u, k, {1, 0, 2});
    //     };
    //     run_mutation_tests(subtest);
    // }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_child_haploid_haploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 2);
    //auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1],potential.shuffle[2]);
    return create_transition_child<1, 1>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_haploid_haploid") {
    auto test = [&](int n, float u, float k, std::vector<int> shuf) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::detail::Potential:: ChildHaploidHaploid, 0, 1, u, 2, 1.1f*u};
        pot.shuffle = shuf;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_haploid_haploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);
        auto mat2 = model.CreateMatrix(n, 1.1f*u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 3);
        REQUIRE(obs.shape(shuffle_pos(0,pot.shuffle)) == n*(n+1)/2);
        REQUIRE(obs.shape(shuffle_pos(1,pot.shuffle)) == n);
        REQUIRE(obs.shape(shuffle_pos(2,pot.shuffle)) == n);
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
                        int o1 = indexes[pot.shuffle[0]];
                        int o2 = indexes[pot.shuffle[1]];
                        int o3 = indexes[pot.shuffle[2]];
                        CHECK(obs(o1,o2,o3) == doctest::Approx(expected));
                    }
                }
            }
        }
    };
    SUBCASE("Shuffle order is {0, 1, 2}.") {
        auto subtest = [&](int n, float u, float k) {
            return test(n, u, k, {0, 1, 2});
        };
        run_mutation_tests(subtest);
    }
    // SUBCASE("Shuffle order is {1, 0, 2}.") {
    //     auto subtest = [&](int n, float u, float k) {
    //         return test(n, u, k, {1, 0, 2});
    //     };
    //     run_mutation_tests(subtest);
    // }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_child_selfing_diploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    //auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);    
    return create_transition_child_selfing<2>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_selfing_diploid") {
    auto test = [&](int n, float u, float k, std::vector<int> shuf) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::detail::Potential::ChildSelfingDiploid, 0, 1, u};
        pot.shuffle = shuf;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_selfing_diploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 2);
        REQUIRE(obs.shape(shuffle_pos(0,pot.shuffle)) == n*(n+1)/2);
        REQUIRE(obs.shape(shuffle_pos(1,pot.shuffle)) == n*(n+1)/2);
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
                        int o1 = indexes[pot.shuffle[0]];
                        int o2 = indexes[pot.shuffle[1]];
                        CHECK(obs(o1,o2) == doctest::Approx(expected));
                    }
                }
            }
        }
    };
    SUBCASE("Shuffle order is {0, 1}.") {
        auto subtest = [&](int n, float u, float k) {
            return test(n, u, k, {0, 1});
        };
        run_mutation_tests(subtest);
    }
    // SUBCASE("Shuffle order is {1, 0}.") {
    //     auto subtest = [&](int n, float u, float k) {
    //         return test(n, u, k, {1, 0});
    //     };
    //     run_mutation_tests(subtest);
    // }
}
// LCOV_EXCL_STOP

template<typename Arg>
mutk::tensor_t create_transition_child_selfing_haploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() == 1);
    //auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);
    return create_transition_child_selfing<1>(model, n, potential, arg);
}

// LCOV_EXCL_START
TEST_CASE("create_transition_child_selfing_haploid") {
    auto test = [&](int n, float u, float k, std::vector<int> shuf) {
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(u);

        potential_t pot{mutk::detail::Potential::ChildSelfingHaploid, 0, 1, u};
        pot.shuffle = shuf;
        KAllelesModel model(k, 0.001, 0, 0, 0);

        auto obs = create_transition_child_selfing_haploid(model, n, pot, mutk::mutation::ANY);
        auto mat1 = model.CreateMatrix(n, u, mutk::mutation::ANY);

        REQUIRE(obs.dimension() == 2);
        REQUIRE(obs.shape(shuffle_pos(0,pot.shuffle)) == n*(n+1)/2);
        REQUIRE(obs.shape(shuffle_pos(1,pot.shuffle)) == n);
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
                    int o1 = indexes[pot.shuffle[0]];
                    int o2 = indexes[pot.shuffle[1]];
                    CHECK(obs(o1,o2) == doctest::Approx(expected));
                }
            }
        }
    };
    SUBCASE("Shuffle order is {0, 1}.") {
        auto subtest = [&](int n, float u, float k) {
            return test(n, u, k, {0, 1});
        };
        run_mutation_tests(subtest);
    }
    // SUBCASE("Shuffle order is {1, 0}.") {
    //     auto subtest = [&](int n, float u, float k) {
    //         return test(n, u, k, {1, 0});
    //     };
    //     run_mutation_tests(subtest);
    // }
}
// LCOV_EXCL_STOP

namespace mutk {
template<typename Arg>
inline
mutk::tensor_t create_mutation_potential(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    using P = mutk::detail::Potential;

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
        return {};
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
