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

#include <mutk/potential.hpp>

using mutk::message_t;

// ==== UNIT POTENTIAL =========================================================

message_t mutk::UnitPotential::Create(size_t n, any_t) {
    auto shape = message_shape(n, labels_);
    return xt::ones<float_t>(shape);
}

message_t mutk::UnitPotential::Create(size_t n, mean_t) {
    auto shape = message_shape(n, labels_);
    return xt::zeros<float_t>(shape);
}

message_t mutk::UnitPotential::Create(size_t n, some_t k) {
    auto shape = message_shape(n, labels_);
    if(+k == 0) {
        return xt::ones<float_t>(shape);
    }
    return xt::zeros<float_t>(shape);
}