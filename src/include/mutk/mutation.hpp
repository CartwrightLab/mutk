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
#include <mutk/potential.hpp>

#include <boost/container/static_vector.hpp>
#include <vector>

namespace mutk {
namespace mutation {

// Tags
struct any_t {};
struct mean_t {};

constexpr any_t ANY = {};
constexpr mean_t MEAN = {};

class Model {
public:
    using tensor_t = mutk::tensor_t;

    virtual tensor_t CreateMatrix(size_t n, float t, any_t) const = 0;
    virtual tensor_t CreateMatrix(size_t n, float t, mean_t) const = 0;
    virtual tensor_t CreateMatrix(size_t n, float t, size_t x) const = 0;

    virtual tensor_t CreatePriorDiploid(size_t n) const = 0;
    virtual tensor_t CreatePriorHaploid(size_t n) const = 0;

    tensor_t CreatePotential(size_t n, const potential_t &potential, any_t);
    tensor_t CreatePotential(size_t n, const potential_t &potential, mean_t);
    tensor_t CreatePotential(size_t n, const potential_t &potential, size_t x);

    virtual ~Model() = default;
};

// A k-alleles model. See Lewis 2001 and Tuffley and Steel 1997.
class KAllelesModel : public Model {
public:
    using tensor_t = Model::tensor_t;

    KAllelesModel(double k, double theta, double hom_bias, double het_bias, double hap_bias) : 
        k_{k}, theta_{theta}, hom_bias_{hom_bias}, het_bias_{het_bias}, hap_bias_{hap_bias} {

        double e = theta/(k-1.0);

        if(k < 2.0) {
            throw std::invalid_argument("Invalid k parameter.");
        }
        if(theta < 0) {
            throw std::invalid_argument("Invalid theta parameter");
        }
        if(hom_bias > 1.0 || hom_bias < -(2.0+e)/((k-1.0)*e)) {
            throw std::invalid_argument("Invalid hom_bias parameter.");
        }
        if(het_bias > 1.0 || het_bias < -(2.0+2.0*e)/((k-2.0)*e)) {
            throw std::invalid_argument("Invalid het_bias parameter.");
        }
        if(hap_bias > 1.0 || hap_bias < -(1.0+e)/((k-1.0)*e)) {
            throw std::invalid_argument("Invalid hap_bias parameter.");
        }
    }

    virtual tensor_t CreateMatrix(size_t n, float t, any_t) const override;
    virtual tensor_t CreateMatrix(size_t n, float t, mean_t) const override;
    virtual tensor_t CreateMatrix(size_t n, float t, size_t x) const override;

    virtual tensor_t CreatePriorDiploid(size_t n) const override;
    virtual tensor_t CreatePriorHaploid(size_t n) const override;

protected:
    double k_;
    double theta_;
    double hom_bias_;
    double het_bias_;
    double hap_bias_;
};

} // namespace mutation
} // namespace mutk

#endif // MUTK_MUTATION_HPP
