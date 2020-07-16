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

// ret(i,j) = P(i|j)
KAllelesModel::matrix_t KAllelesModel::CreateMatrix(int n, float t, any_t) const {
    assert(n > 0);
    assert(n <= 5);

    double beta = t*k_/(k_-1.0);
    double p_ji = -1.0/k_*expm1(-beta);
    double p_jj = exp(-beta) + p_ji;

    Tensor<2> ret(n,n);

    return ret.generate([&](const auto& coords) {
        int i = coords[0];
        int j = coords[1];
        return (i == j) ? p_jj : p_ji;
    });
}

TEST_CASE("[libmutk] mutation::KAllelesModel::CreateMatrix with any_t") {
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
        REQUIRE(obs.dimensions().size() == 2);
        REQUIRE(obs.dimension(0) == n);
        REQUIRE(obs.dimension(1) == n);
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

// ret(i,j) = P(i & x mutations | j)
KAllelesModel::matrix_t KAllelesModel::CreateMatrix(int n, float t, int x) const {
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

    Tensor<2> ret(n,n);

    return ret.generate([&](const auto& coords) {
        int i = coords[0];
        int j = coords[1];
        return (i == j) ? p_x*p_jj : p_x*p_ji;
    });
}

TEST_CASE("[libmutk] mutation::KAllelesModel::CreateMatrix with int") {
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
            REQUIRE(obs.dimensions().size() == 2);
            REQUIRE(obs.dimension(0) == n);
            REQUIRE(obs.dimension(1) == n);
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


// ret(i,j) = E[num of mutations | i,j]*P(i|j)
KAllelesModel::matrix_t KAllelesModel::CreateMatrix(int n, float t, mean_t) const {
    assert(n > 0);
    assert(n <= 5);

    Tensor<2> ret(n,n);

    double beta = k_*t/(k_-1.0);
    double p_jj = -t/k_*expm1(-beta);
    double p_ji = (t-p_jj)/(k_-1.0);

    return ret.generate([&](const auto& coords) {
        int i = coords[0];
        int j = coords[1];
        return (i == j) ? p_jj : p_ji;
    });

    return ret;
}

TEST_CASE("[libmutk] mutation::KAllelesModel::CreateMatrix with mean_t") {
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
        REQUIRE(obs.dimensions().size() == 2);
        REQUIRE(obs.dimension(0) == n);
        REQUIRE(obs.dimension(1) == n);
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

KAllelesModel::tensor_t KAllelesModel::CreatePriorHaploid(int n) const {
    double k = k_;
    double e = theta_/(k-1.0);

    double p_R = (1.0+e+(k-1.0)*e*hap_bias_)/(1.0+k*e);
    double p_A = (e-e*hap_bias_)/(1.0+k*e);

    tensor_t ret(mutk::dim_width<1>(n));
    return ret.generate([&](const auto &coords) {
        Eigen::DenseIndex i = coords[0];
        return (i == 0) ? p_R : p_A;
    });
}

TEST_CASE("[libmutk] mutation::KAllelesModel::CreatePriorHaploid") {
    auto test_haploid = [](int n, float theta, float hap_bias,
        float k) {
        CAPTURE(n);
        CAPTURE(theta);
        CAPTURE(hap_bias);
        CAPTURE(k);

        KAllelesModel model(k, theta, 0, 0, hap_bias);

        auto obs = model.CreatePriorHaploid(n);

        if(n == k) {
            mutk::Tensor<0> s = obs.sum();
            CHECK(s(0) == doctest::Approx(1.0f));
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

KAllelesModel::tensor_t KAllelesModel::CreatePriorDiploid(int n) const {
    double k = k_;
    double e = theta_/(k-1.0);

    double p_hom = (1.0+e)/(1.0+k*e);
    double p_hetk = e/(1.0+k*e);

    double p_RR = p_hom*(2.0+e+(k-1.0)*e*hom_bias_)/(2.0+k*e);
    double p_AA = p_hom*(e-e*hom_bias_)/(2.0+k*e);

    double p_RA = p_hetk*(2.0+2.0*e+(k-2.0)*e*het_bias_)/(2.0+k*e);
    double p_AB = p_hetk*(2.0*e-2.0*e*het_bias_)/(2.0+k*e);

    tensor_t ret(mutk::dim_width<2>(n));
    return ret.generate([&](const auto &coords){
        Eigen::DenseIndex a = ALLELE[coords[0]][0];
        Eigen::DenseIndex b = ALLELE[coords[0]][1];
        if(a == b) {
            return (a == 0) ? p_RR : p_AA;
        }
        return (a == 0) ? p_RA : p_AB;
    });
}

TEST_CASE("[libmutk] mutation::KAllelesModel::CreatePriorDiploid") {
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
            mutk::Tensor<0> s = obs.sum();
            CHECK(s(0) == doctest::Approx(1.0f));
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

template<typename Arg>
mutk::Tensor<2> create_transition_clone_haploid_impl(const mutk::mutation::Model &model,
    int n, float t, Arg arg) {

    return model.CreateMatrix(n, t, arg);
}

template<typename Arg>
mutk::Tensor<2> create_transition_gamete_diploid_impl(const mutk::mutation::Model &model, 
    int n, float t, Arg arg) {

    auto mat = create_transition_clone_haploid_impl(model, n, t, arg);

    mutk::Tensor<2> ret(mutk::dim_width<1>(n), mutk::dim_width<2>(n));

    // x/y -> a
    return ret.generate([&](const auto& coords) {
        int a = coords[0];
        int x = ALLELE[coords[1]][0];
        int y = ALLELE[coords[1]][1];
        return 0.5*(mat(a,x)+mat(a,y));
    });    
}

template<int N, typename Arg>
inline
mutk::Tensor<2> create_transition_haploid(const mutk::mutation::Model &model, 
    int n, float t, Arg arg) {
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

    using array_t = Eigen::array<Eigen::DenseIndex, TENSOR_SIZE>;

    static array_t Dims(int n) {
        return {mutk::dim_width<2>(n), mutk::dim_width<2>(n)};
    }
    static Eigen::DenseIndex G1(const array_t &coords) {
        return ALLELE[coords[1]][0];
    }
    static Eigen::DenseIndex G2(const array_t &coords) {
        return ALLELE[coords[1]][1];
    }
};

// G1 x G2 => a/b
template<int N, int M>
struct transition_traits<child_diploid_tag<N, M>> {
    static constexpr int MAT1_SIZE = N;
    static constexpr int MAT2_SIZE = M;
    static constexpr int TENSOR_SIZE = 3;

    using array_t = Eigen::array<Eigen::DenseIndex, TENSOR_SIZE>;

    static array_t Dims(int n) {
        return {mutk::dim_width<2>(n), mutk::dim_width<N>(n), mutk::dim_width<M>(n)};
    }
    static Eigen::DenseIndex G1(const array_t &coords) {
        return coords[1];
    }
    static Eigen::DenseIndex G2(const array_t &coords) {
        return coords[2];
    }
};

// G1 x G1 => a/b
template<int N>
struct transition_traits<child_selfing_tag<N>> {
    static constexpr int MAT1_SIZE = N;
    static constexpr int MAT2_SIZE = N;
    static constexpr int TENSOR_SIZE = 2;

    using array_t = Eigen::array<Eigen::DenseIndex, TENSOR_SIZE>;

    static array_t Dims(int n) {
        return {mutk::dim_width<2>(n), mutk::dim_width<N>(n)};
    }
    static Eigen::DenseIndex G1(const array_t &coords) {
        return coords[1];
    }
    static Eigen::DenseIndex G2(const array_t &coords) {
        return coords[1];
    }
};

template<typename Tag>
struct generate_transition_diploid {
    using traits = transition_traits<Tag>;

    template<typename Arg1, typename Arg2>
    inline
    mutk::Tensor<traits::TENSOR_SIZE>
    operator()(const mutk::mutation::Model &model, 
        int n, float t1, float t2, Arg1 arg1, Arg2 arg2) const {
     
        auto mat1 = create_transition_haploid<traits::MAT1_SIZE>(model, n, t1, arg1);
        auto mat2 = create_transition_haploid<traits::MAT2_SIZE>(model, n, t2, arg2);

        mutk::Tensor<traits::TENSOR_SIZE> ret(traits::Dims(n));

        // G1/G2 => a/b OR G1 x G2 => a/b
        return ret.generate([&](const auto& coords) {
            Eigen::DenseIndex a = ALLELE[coords[0]][0];
            Eigen::DenseIndex b = ALLELE[coords[0]][1];
            Eigen::DenseIndex g1 = traits::G1(coords);
            Eigen::DenseIndex g2 = traits::G2(coords);

            float f = mat1(a,g1) * mat2(b,g2);
            if(a != b) {
                f += mat1(b,g1) * mat2(a,g2);
            }
            return f;            
        });
    }
};

template<typename Tag, typename Arg>
auto create_transition_diploid_impl(const mutk::mutation::Model &model, 
    int n, float t1, float t2, Arg arg) {
    generate_transition_diploid<Tag> gen;

    if constexpr(std::is_same_v<Arg, mutk::mutation::any_t>) {
        return gen(model, n, t1, t2, arg, arg);
    } else if constexpr(std::is_same_v<Arg, mutk::mutation::mean_t>) {
        return gen(model, n, t1, t2, arg, mutk::mutation::ANY)
                + gen(model, n, t1, t2, mutk::mutation::ANY, arg);
    } else if constexpr(std::is_integral_v<Arg>) {
        auto ret = gen(model, n, t1, t2, 0, arg);
        for(Arg a=1; a < arg; ++a) {
            ret += gen(model, n, t1, t2, a, arg-a);
        }
        return ret;
    } else {
        static_assert(std::is_integral_v<Arg>
            || std::is_same_v<Arg, mutk::mutation::mean_t>
            || std::is_same_v<Arg, mutk::mutation::any_t> );
    }     
}


template<typename Arg>
mutk::Tensor<2> create_transition_clone_diploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() >= 1);
    auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);
    return create_transition_diploid_impl<clone_diploid_tag>(model, n, potential.parents[0].second,
        potential.parents[0].second, arg).shuffle(shuffle_axes);
}

template<int N, int M, typename Arg>
inline
mutk::Tensor<3> create_transition_child(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() >= 2);
    auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1],potential.shuffle[2]);
    return create_transition_diploid_impl<child_diploid_tag<N,M>>(model, n, potential.parents[0].second,
        potential.parents[1].second, arg).shuffle(shuffle_axes);
}

template<int N, typename Arg>
inline
mutk::Tensor<2> create_transition_child_selfing(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() >= 1);
    auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);
    return create_transition_diploid_impl<child_selfing_tag<N>>(model, n, potential.parents[0].second,
        potential.parents[0].second, arg).shuffle(shuffle_axes);
}

template<typename Arg>
mutk::Tensor<2> create_transition_clone_haploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() >= 1);
    auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);
    return create_transition_haploid<1>(model, n, potential.parents[0].second, arg).shuffle(shuffle_axes);
}

template<typename Arg>
mutk::Tensor<2> create_transition_gamete_diploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() >= 1);
    auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);
    return create_transition_haploid<2>(model, n, potential.parents[0].second, arg).shuffle(shuffle_axes);
}

template<typename Arg>
mutk::Tensor<3> create_transition_child_diploid_diploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() >= 2);
    auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1],potential.shuffle[2]);
    return create_transition_child<2, 2>(model, n, potential, arg).shuffle(shuffle_axes);
}

template<typename Arg>
mutk::Tensor<3> create_transition_child_haploid_diploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() >= 2);
    auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1],potential.shuffle[2]);
    return create_transition_child<1, 2>(model, n, potential, arg).shuffle(shuffle_axes);
}

template<typename Arg>
mutk::Tensor<3> create_transition_child_diploid_haploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() >= 2);
    auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1],potential.shuffle[2]);
    return create_transition_child<2, 1>(model, n, potential, arg).shuffle(shuffle_axes);
}

template<typename Arg>
mutk::Tensor<3> create_transition_child_haploid_haploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() >= 2);
    auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1],potential.shuffle[2]);
    return create_transition_child<1, 1>(model, n, potential, arg).shuffle(shuffle_axes);
}

template<typename Arg>
mutk::Tensor<2> create_transition_child_selfing_diploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() >= 1);
    auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);    
    return create_transition_child_selfing<2>(model, n, potential, arg).shuffle(shuffle_axes);
}

template<typename Arg>
mutk::Tensor<2> create_transition_child_selfing_haploid(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    assert(potential.parents.size() >= 1);
    auto shuffle_axes = mutk::tensor_dims(potential.shuffle[0],potential.shuffle[1]);
    return create_transition_child_selfing<1>(model, n, potential, arg).shuffle(shuffle_axes);
}

namespace mutk {
template<typename Arg>
inline
mutk::Tensor<1> create_mutation_potential(const mutk::mutation::Model &model, 
    int n, const potential_t &potential, Arg arg) {
    using P = mutk::detail::Potential;

    int g = dim_width<2>(n);

    switch(potential.type) {
    case P::CloneDiploid:
        return create_transition_clone_diploid(model, n, potential, arg)
            .reshape(tensor_dims(g*g));
    case P::CloneHaploid:
        return create_transition_clone_haploid(model, n, potential, arg)
            .reshape(tensor_dims(n*n));
    case P::GameteDiploid:
        return create_transition_gamete_diploid(model, n, potential, arg)
            .reshape(tensor_dims(n*g));
    case P::ChildDiploidDiploid:
        return create_transition_child_diploid_diploid(model, n, potential, arg)
            .reshape(tensor_dims(g*g*g));
    case P::ChildHaploidDiploid:
        return create_transition_child_haploid_diploid(model, n, potential, arg)
            .reshape(tensor_dims(g*n*g));
    case P::ChildDiploidHaploid:
        return create_transition_child_diploid_haploid(model, n, potential, arg)
            .reshape(tensor_dims(g*g*n));
    case P::ChildHaploidHaploid:
        return create_transition_child_haploid_haploid(model, n, potential, arg)
            .reshape(tensor_dims(g*n*n));
    case P::ChildSelfingDiploid:
        return create_transition_child_selfing_diploid(model, n, potential, arg)
            .reshape(tensor_dims(g*g));
    case P::ChildSelfingHaploid:
        return create_transition_child_selfing_haploid(model, n, potential, arg)
            .reshape(tensor_dims(g*n));
    case P::Unit:
        return {};
    default:
        break;
    };
    return {};
}
}

mutk::Tensor<1> mutk::mutation::Model::CreatePotential(int n, const potential_t &potential, any_t arg) {
    return create_mutation_potential(*this, n, potential, arg);
}

mutk::Tensor<1> mutk::mutation::Model::CreatePotential(int n, const potential_t &potential, mean_t arg) {
    return create_mutation_potential(*this, n, potential, arg);
}

mutk::Tensor<1> mutk::mutation::Model::CreatePotential(int n, const potential_t &potential, int arg) {
    return create_mutation_potential(*this, n, potential, arg);
}
