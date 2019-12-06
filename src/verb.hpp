/*
# Copyright (c) 2019 Reed A. Cartwright <reed@cartwright.ht>
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

#ifndef MUTK_VERB_HPP
#define MUTK_VERB_HPP

#include <iostream>

#include <mutk/mutk.hpp>

namespace mutk {
namespace verb {

inline
int check_version_number() {
    if(mutk::version_number_check_equal(MUTK_VERSION_INTEGER) == false) {
        std::cerr << "ERROR: Version mismatch between headers (#"
                  << MUTK_VERSION_INTEGER << ") and library (#"
                  << mutk::version_integer() << ").\n";
        std::cerr << "       MUTK linked against wrong version of library.\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

#define MUTK_VERB_RUNTIME_CHECK_VERSION_NUMBER_OR_RETURN() \
    do { \
        auto check = mutk::verb::check_version_number(); \
        if(check != EXIT_SUCCESS) { \
            return check; \
        } \
    } while(false) \

} // namespace verb
} // namespace mutk


#endif