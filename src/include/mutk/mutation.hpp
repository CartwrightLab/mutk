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

#ifndef MUTK_MUTATION_HPP
#define MUTK_MUTATION_HPP

#include <stdexcept>

#include "message.hpp"

namespace mutk {

// A k-alleles model. See Lewis 2001 and Tuffley and Steel 1997.
//
// P(j|i) = 1/k * (1-exp(-beta*t)) + I(j == i)*exp(-beta*t)
//    where beta = k/(k-1)
//    expected number of mutations = t
//    expected number of events = beta*t
//
//    1/k mutations are virtual and (k-1)/k are real
//
//    Given at least 1 event, the expected number of events is 
//      (beta*t)/(1-exp(-beta*t))
//
//    And expected number of mutations is t/(1-exp(-beta*t))
//
class MutationModel {
public:
    using array_t = mutk::message_t;

    MutationModel(float_t k, float_t theta, float_t hom_bias, float_t het_bias, float_t hap_bias) : 
        k_{k}, theta_{theta}, hom_bias_{hom_bias}, het_bias_{het_bias}, hap_bias_{hap_bias} {

        float_t e = theta/(k-1.0f);

        if(k < 2.0f) {
            throw std::invalid_argument("Invalid k parameter.");
        }
        if(theta < 0.0f) {
            throw std::invalid_argument("Invalid theta parameter");
        }
        if(hom_bias > 1.0f || hom_bias < -(2.0f+e)/((k-1.0f)*e)) {
            throw std::invalid_argument("Invalid hom_bias parameter.");
        }
        if(het_bias > 1.0f || het_bias < -(2.0f+2.0f*e)/((k-2.0f)*e)) {
            throw std::invalid_argument("Invalid het_bias parameter.");
        }
        if(hap_bias > 1.0f || hap_bias < -(1.0f+e)/((k-1.0f)*e)) {
            throw std::invalid_argument("Invalid hap_bias parameter.");
        }
    }

    array_t CreateTransitionMatrix(message_size_t n, float_t t) const;
    array_t CreateMeanMatrix(message_size_t n, float_t t) const;
    array_t CreateCountMatrix(message_size_t n, float_t t, int x) const;

protected:
    float_t k_;
    float_t theta_;
    float_t hom_bias_;
    float_t het_bias_;
    float_t hap_bias_;
};

} // namespace mutk

#endif // MUTK_MUTATION_HPP
