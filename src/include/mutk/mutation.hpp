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

#ifndef MUTK_MUTATION_HPP
#define MUTK_MUTATION_HPP

#include <mutk/memory.hpp>

namespace mutk {
namespace mutation {

class Model {
public:
    using tensor_t = mutk::Tensor<2>;
    virtual tensor_t TransitionMatrix(float t, int n) const = 0;
    virtual tensor_t EventTransitionMatrix(float t, int n, int x) const = 0;
    virtual tensor_t MeanTransitionMatrix(float t, int n) const = 0;
};

// A k-alleles model. See Lewis 2001 and Tuffley and Steel 1997.
class KAllelesModel : public Model {
public:
    using tensor_t = Model::tensor_t;

    KAllelesModel(double k) : k_{k} {
        assert(k_ >= 2.0);        
    }

    virtual tensor_t TransitionMatrix(float t, int n) const override;
    virtual tensor_t EventTransitionMatrix(float t, int n, int x) const override;
    virtual tensor_t MeanTransitionMatrix(float t, int n) const override;

protected:
    double k_;
};



/*
struct transition_t {};
struct mean_t {};

inline
tensor_t mitosis_haploid_matrix(int size, Model m, transition_t) {
    return m.TransitionMatrix(size);
}

inline
tensor_t mitosis_haploid_matrix(int size, Model m, mean_t) {
    return m.MeanTransitionMatrix(size);
}

inline
tensor_t mitosis_haploid_matrix(int size, Model m, int count) {
    return m.EventTransitionMatrix(size, count);
}

namespace detail {
inline
void mitosis_diploid_matrix_op(const tensor_t& matA, const tensor_t& matB, tensor_t *p) {
    assert(matA.cols() == matA.rows() && matA.cols() > 0);
    assert(matB.cols() == matB.rows() && matB.cols() > 0);
    assert(matA.cols() == matB.cols());

    const int num_alleles = matB.cols();
    const int num_genotypes = num_alleles*(num_alleles+1)/2;
    assert(p != nullptr && p->cols() == p->rows() && p->cols() == num_genotypes);

    for(int a=0,i=0; a< num_alleles; ++a) { // loop over all child genotypes 
        for(int b=0; b<=a; ++b, ++i) { 
            for(int x=0,j=0; x<num_alleles; ++x) { // loop over all parent genotypes
                for(int y=0; y<=x; ++y,++j) {
                    // x/y => a/b
                    // combine the results from the two phases unless a == b
                    (*p)(j,i) += matA(x,a) * matB(y,b);
                    if(a != b) {
                        (*p)(j,i) += matA(x,b) * matB(y,a);
                    }
                }
            }
        }
    }
}

} // namespace detail

inline
tensor_t mitosis_diploid_matrix(int size, Model m, transition_t) {
    assert(size > 0);
    const int num_genotypes = size*(size+1)/2;
    
    tensor_t ret = tensor_t::Zero(num_genotypes, num_genotypes);

    auto mat = mitosis_haploid_matrix(size, m, transition_t{});
    detail::mitosis_diploid_matrix_op(mat,mat,&ret);

    return ret;
}

inline
tensor_t mitosis_diploid_matrix(int size, Model m, int count) {
    assert(size > 0);
    const int num_genotypes = size*(size+1)/2;

    tensor_t ret = tensor_t::Zero(num_genotypes, num_genotypes);
    
    for(int n=0; n<=count; ++n) {
        auto mat1 = mitosis_haploid_matrix(size, m, n);
        auto mat2 = mitosis_haploid_matrix(size, m, count-n);
        detail::mitosis_diploid_matrix_op(mat1,mat2,&ret);
    }

    return ret;
}

inline
tensor_t mitosis_diploid_matrix(int size, Model m, mean_t) {
    assert(size > 0);
    const int num_alleles = size;
    const int num_genotypes = num_alleles*(num_alleles+1)/2;
    
    tensor_t ret = tensor_t::Zero(num_genotypes, num_genotypes);

    auto mat = mitosis_haploid_matrix(size, m, transition_t{});
    auto avg = mitosis_haploid_matrix(size, m, mean_t{});

    detail::mitosis_diploid_matrix_op(mat, avg, &ret);
    detail::mitosis_diploid_matrix_op(avg, mat, &ret);

    return ret;
}

namespace detail {
inline
void meiosis_haploid_matrix_op(const tensor_t& matA, tensor_t *p) {
    assert(matA.cols() == matA.rows() && matA.cols() > 0);

    const int num_alleles = matA.cols();
    const int num_genotypes = num_alleles*(num_alleles+1)/2;
    assert(p != nullptr && p->cols() == num_alleles && p->rows() == num_genotypes);

    for(int i=0; i<num_alleles; ++i) { // loop over all child haplotypes
        for(int x=0,j=0; x<num_alleles; ++x) { // loop over all parent genotypes
            for(int y=0; y<=x; ++y,++j) {
                (*p)(j,i) += 0.5*(matA(x,i)+matA(y,i));
            }
        }
    }
}
} // namespace detail

template<typename T>
inline
tensor_t meiosis_haploid_matrix(int size, Model m, T arg) {
    assert(size > 0);
    const int num_genotypes = size*(size+1)/2;
    
    tensor_t ret = tensor_t::Zero(num_genotypes, size);

    auto mat = mitosis_haploid_matrix(size, m, arg);
    detail::meiosis_haploid_matrix_op(mat, &ret);

    return ret;
}

template<typename T>
inline
tensor_t mitosis_matrix(int size, Model m, T arg, int parent_ploidy) {
    assert(parent_ploidy == 1 || parent_ploidy == 2);
    if(parent_ploidy == 1) {
        return mitosis_haploid_matrix(size, m, arg);
    }
    return mitosis_diploid_matrix(size, m, arg);    
}

template<typename T>
inline
tensor_t gamete_matrix(int size, Model m, T arg, int parent_ploidy) {
    assert(parent_ploidy == 1 || parent_ploidy == 2);
    if(parent_ploidy == 1) {
        return mitosis_haploid_matrix(size, m, arg);
    }
    return meiosis_haploid_matrix(size, m, arg);
}

inline
int number_of_parent_genotypes(const int num_alleles, const int ploidy) {
    assert(ploidy == 1 || ploidy == 2);
    return ((ploidy == 1) ? num_alleles : num_alleles*(num_alleles+1)/2);
}

inline
int number_of_parent_genotype_pairs(const int num_alleles, const int dad_ploidy, const int mom_ploidy) {
    return number_of_parent_genotypes(num_alleles, dad_ploidy)
        * number_of_parent_genotypes(num_alleles, mom_ploidy);
}

namespace detail {
inline
void meiosis_matrix_op(const tensor_t& matA, const tensor_t& matB, tensor_t *p) {
    assert(matA.cols() == matB.cols());

    const int num_alleles = matA.cols();
    const int num_genotypes = num_alleles*(num_alleles+1)/2;
    assert(p != nullptr);
    assert(p->cols() == num_genotypes);
    assert(p->rows() == matA.rows()*matB.rows() );

    for(int a=0,i=0; a < matA.cols(); ++a) { // loop over all child genotypes 
        for(int b=0; b <= a; ++b, ++i) {
            for(int x=0,j=0; x < matA.rows(); ++x) { // loop over all parent genotypes
                for(int y=0; y < matB.rows(); ++y,++j) {
                    // fold the tensor_t
                    (*p)(j,i) += matA(x,a) * matB(y,b);
                    if(b != a) {
                        (*p)(j,i) += matA(x,b) * matB(y,a);
                    }
                }
            }
        }
    }
}
} // detail

inline
tensor_t meiosis_matrix(int size, Model dad_m, Model mom_m, transition_t, int dad_ploidy, int mom_ploidy) {
    assert(dad_ploidy == 1 || dad_ploidy == 2);
    assert(mom_ploidy == 1 || mom_ploidy == 2);

    const int num_alleles = size;
    const int num_genotypes = num_alleles*(num_alleles+1)/2;

    // Construct Mutation Process
    tensor_t ret = tensor_t::Zero(number_of_parent_genotype_pairs(num_alleles,
            dad_ploidy, mom_ploidy), num_genotypes);
    auto dad = gamete_matrix(size, dad_m, transition_t{}, dad_ploidy);
    auto mom = gamete_matrix(size, mom_m, transition_t{}, mom_ploidy);

    detail::meiosis_matrix_op(dad,mom,&ret);
    return ret;
}

inline
tensor_t meiosis_matrix(int size, Model dad_m, Model mom_m, mean_t, int dad_ploidy, int mom_ploidy) {
    assert(dad_ploidy == 1 || dad_ploidy == 2);
    assert(mom_ploidy == 1 || mom_ploidy == 2);

    const int num_alleles = size;
    const int num_genotypes = num_alleles*(num_alleles+1)/2;

    // Construct Mutation Process
    tensor_t ret = tensor_t::Zero(number_of_parent_genotype_pairs(num_alleles,
            dad_ploidy, mom_ploidy), num_genotypes);

    auto dad = gamete_matrix(size, dad_m, transition_t{}, dad_ploidy);
    auto dad_mean = gamete_matrix(size, dad_m, mean_t{}, dad_ploidy);
    auto mom = gamete_matrix(size, mom_m, transition_t{}, mom_ploidy);
    auto mom_mean = gamete_matrix(size, mom_m, mean_t{}, mom_ploidy);

    detail::meiosis_matrix_op(dad,mom_mean,&ret);
    detail::meiosis_matrix_op(dad_mean,mom,&ret);    
    return ret;
}

inline
tensor_t meiosis_matrix(int size, Model dad_m, Model mom_m, int count, int dad_ploidy, int mom_ploidy) {
    assert(dad_ploidy == 1 || dad_ploidy == 2);
    assert(mom_ploidy == 1 || mom_ploidy == 2);

    const int num_alleles = size;
    const int num_genotypes = num_alleles*(num_alleles+1)/2;

    // Construct Mutation Process
    tensor_t ret = tensor_t::Zero(number_of_parent_genotype_pairs(num_alleles,
            dad_ploidy, mom_ploidy), num_genotypes);

   for(int n=0; n<=count; ++n) {
        auto dad = gamete_matrix(size, dad_m, n, dad_ploidy);
        auto mom = gamete_matrix(size, mom_m, count-n, mom_ploidy);
     
        detail::meiosis_matrix_op(dad,mom,&ret);
    }
    return ret;
}

// k-alleles model from Watterson and Guess (1977) https://doi.org/10.1016/0040-5809(77)90023-5

inline
tensor_t population_prior_diploid(int num_obs_alleles, double theta, double hom_bias, double het_bias,
    double kalleles) {
    assert(num_obs_alleles >= 0);

    double k = kalleles;
    double e = theta/(k-1.0);

    double p_hom = (1.0+e)/(1.0+k*e);
    double p_hetk = e/(1.0+k*e);

    double p_RR = p_hom*(2.0+e+(k-1.0)*e*hom_bias)/(2.0+k*e);
    double p_AA = p_hom*(e-e*hom_bias)/(2.0+k*e);

    double p_RA = p_hetk*(2.0+2.0*e+(k-2.0)*e*het_bias)/(2.0+k*e);
    double p_AB = p_hetk*(2.0*e-2.0*e*het_bias)/(2.0+k*e);

    tensor_t ret{num_obs_alleles*(num_obs_alleles+1)/2};

    int n=0;
    for(int i=0;i<num_obs_alleles;++i) {
        for(int j=0;j<i;++j) {
            ret(n++) = (j==0 || i==0) ? p_RA : p_AB;
        }
        ret(n++) = (i==0) ? p_RR : p_AA;
    }

    return ret;
}

inline
tensor_t population_prior_haploid(int num_obs_alleles, double theta, double hap_bias,
    double kalleles) {
    assert(num_obs_alleles >= 1);

    double k = kalleles;
    double e = theta/(k-1.0);

    double p_R = (1.0+e+(k-1.0)*e*hap_bias)/(1.0+k*e);
    double p_A = (e-e*hap_bias)/(1.0+k*e);

    tensor_t ret{num_obs_alleles};
    ret(0) = p_R;
    for(int n=1;n<num_obs_alleles;++n) {
        ret(n) = p_A;
    }

    return ret;
}

inline bool population_prior_check(double theta, double hom_bias, double het_bias, double hap_bias, double kalleles) {
    double k = kalleles;
    double e = theta/(k-1.0);

    if(e < 0) {
        return false;
    }
    if(k < 2.0) {
        return false;
    }
    if(hom_bias > 1.0 || hom_bias < -(2.0+e)/((k-1.0)*e)) {
        return false;
    }
    if(het_bias > 1.0 || het_bias < -(2.0+2.0*e)/((k-2.0)*e)) {
        return false;
    }
    if(hap_bias > 1.0 || hap_bias < -(1.0+e)/((k-1.0)*e)) {
        return false;
    }
    return true;
}

*/

} // namespace mutation
} // namespace mutk

#endif // MUTK_MUTATION_HPP
