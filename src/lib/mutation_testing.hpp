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

#ifndef MUTK_MUTATION_TESTING_HPP
#define MUTK_MUTATION_TESTING_HPP

#include <mutk/message.hpp>

// Libraries needed for testing
#include <boost/numeric/ublas/matrix.hpp>

// Structure useful for unit testing
struct kalleles_test_mat {
    using mat_t = boost::numeric::ublas::matrix<float,
        boost::numeric::ublas::column_major,std::vector<float>>;

    mat_t Q;

    kalleles_test_mat(size_t n, float k) : Q(n+1,n+1) {
        using namespace boost::numeric::ublas;
        // Build Q matrix from Tuffley and Steel (1997)
        vector<float> freqs(n+1, 1.0);
        freqs[n] = k-n;

        for(size_t i=0; i<=n; ++i) {
            for(size_t j=0; j<=n; ++j) {
                Q(i,j) = freqs[j];
            }
            Q(i,i) = Q(i,i)-k;
        }

        // Scale matrix into substitution time
        Q *= 1.0/(k-1.0);
    }
};

template<typename F>
void run_mutation_tests(F test) {
    test(2, 0.0,  0.0, 4.0);
    test(4, 0.0,  1e-8, 4.0);

    test(1, 1e-8, 4e-8, 4.0);
    test(2, 1e-8, 3e-8, 4.0);
    test(3, 1e-8, 2e-8, 4.0);
    test(4, 1e-8, 1e-8, 4.0);
    test(4, 1e-3, 1e-6, 4.0);
    test(4, 1e-6, 1e-3, 4.0);

    test(4, 1e-8, 1e-6, 4.0);
    test(4, 1e-9, 1e-6, 5.0);
    test(4, 1e-6, 1e-6, 6.0);
    test(4, 1e-5, 1e-6, 7.0);
}

// for debugging purposes
namespace std {
inline
std::ostream& operator<< (std::ostream& os, const mutk::message_t::shape_type & value) {
    os << "{";
    for(size_t i = 0; i < value.size(); ++i) {
        if(i > 0) {
            os << ", ";
        }
        os << value[i];
    }
    os << "}";
    return os;
}
}

#endif
