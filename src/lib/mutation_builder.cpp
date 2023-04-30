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

#include <cmath>
#include <iterator>
#include <utility>
#include <mutk/message.hpp>
#include <mutk/mutation.hpp>

#include <boost/math/special_functions/binomial.hpp>

TEST_CASE("MutationMessageBuilder") {
    using Semiring = mutk::mutation_semiring::Probability;
    using Builder = mutk::MutationMessageBuilder<Semiring>;

    using S = mutk::message_t::shape_type; S{};

    {
        Builder builder({2,2});

        builder.AddTransition(0, 0, 1.0, Semiring(4,0.001));
        builder.AddTransition(1, 1, 1.0, Semiring(4,0.001));

        auto msg = builder.Create(4);

        double k = 4;
        double u = 0.001;
        double beta = k/(k-1.0);
        double pij = -1.0/k*expm1(-beta*u);
        double pii = exp(-beta*u) + pij;

        auto expected = mutk::message_t::from_shape({10,10});

        for(int a=0,i=0;a<4;++a) {
            for(int b=0;b<=a;++b,++i) {
                for(int x=0,j=0;x<4;++x) {
                    for(int y=0;y<=x;++y,++j) {
                        // ab -> xy
                        double ax = (a == x) ? pii : pij;
                        double by = (b == y) ? pii : pij;
                        double ay = (a == y) ? pii : pij;
                        double bx = (b == x) ? pii : pij;
                        double ret = ax*by;
                        if(x != y) {
                            ret += ay*bx;
                        }
                        expected(i,j) = ret;
                    }
                }
            }
        }

        CHECK(msg.shape() == expected.shape());
        CHECK_EQ_RANGES(msg, expected);
    }
    {
        Builder builder({2,2,0});

        auto expected = mutk::message_t::from_shape({10,10});
        expected.fill(1.0);
        auto msg = builder.Create(4);
        CHECK(msg.shape() == expected.shape());
        CHECK_EQ_RANGES(msg, expected);
    }
    {
        Builder builder({2,1,2});

        builder.AddTransition(0, 0, 0.5, Semiring(4,0.0001));
        builder.AddTransition(0, 1, 0.5, Semiring(4,0.0001));
        builder.AddTransition(1, 2, 1.0, Semiring(4,0.001));

        auto msg = builder.Create(2);

        double k = 4;
        double u = 0.0001;
        double v = 0.001;
        double beta = k/(k-1.0);
        double piju = -1.0/k*expm1(-beta*u);
        double piiu = exp(-beta*u) + piju;
        double pijv = -1.0/k*expm1(-beta*v);
        double piiv = exp(-beta*v) + pijv;

        auto expected = mutk::message_t::from_shape({3,2,3});

        int n = 2;
        for(int a=0,i=0;a<n;++a) {
            for(int b=0;b<=a;++b,++i) {
                for(int c=0,j=0;c<n;++c,++j) {
                    for(int x=0,k=0;x<n;++x) {
                        for(int y=0;y<=x;++y,++k) {
                            // ab -> x * c -> y
                            double ax = (a == x) ? piiu : piju;
                            double bx = (b == x) ? piiu : piju;
                            double abx = 0.5*(ax+bx);
                            double cy = (c == y) ? piiv : pijv;
                            double ret = abx*cy;
                            if(x != y) {
                                ax = (a == y) ? piiu : piju;
                                bx = (b == y) ? piiu : piju;
                                abx = 0.5*(ax+bx);
                                cy = (c == x) ? piiv : pijv;
                                ret += abx*cy;
                            }
                            expected(i,j,k) = ret;
                        }
                    }
                }
            }
        }

        // 00 x 0 -> 01

        CHECK(msg.shape() == expected.shape());
        CHECK_EQ_RANGES(msg, expected);
    }

    //std::cout << msg << std::endl;
}
