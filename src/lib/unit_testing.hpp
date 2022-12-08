/*
# Copyright (c) 2021 Reed A. Cartwright <racartwright@gmail.com>
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
#include <sstream>

namespace fragmites {
template<typename A, typename B>
inline
void check_eq_ranges_impl(const A& a, const B& b, const char* msg, const char* file, int line) {
    std::ostringstream oss;
    if(std::size(a) == std::size(b)) {
        auto ait = std::begin(a);
        auto bit = std::begin(b);
        size_t bad = 0;
        for(size_t i = 0; i < std::size(a); ++i, ++ait, ++bit) {
            if(*ait == *bit)
                continue;
            if(++bad < 5) {
                oss << "\n  pos " << i << ": " << *ait << " != " << *bit;
            }
        }
        if(bad == 0) {
            return;
        }
        if(bad >= 5) {
            oss << "\n  " << bad << " positions in total NOT correct.";
        }
    } else {
        oss << "\n  size: " << std::size(a) << " != " << std::size(b);
    }
    ADD_FAIL_CHECK_AT(file, line, doctest::Color::Cyan <<
        "CHECK_EQ_RANGES( " << std::string(msg) << " )" << doctest::Color::None << " is NOT correct!"
        << oss.str());
}
}

#define CHECK_EQ_RANGES(x, y) fragmites::check_eq_ranges_impl(x, y, #x ", " #y, __FILE__, __LINE__)
