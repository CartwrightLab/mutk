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

#include <doctest/doctest.h>

#include <mutk/mutk.hpp>

// do some version number sanity checks
static_assert(MUTK_VERSION_MAJOR >= 0 && MUTK_VERSION_MAJOR < 1000,  // NOLINT
              "MUTK major version must be less than 1000.");
static_assert(MUTK_VERSION_MINOR >= 0 && MUTK_VERSION_MINOR < 1000,  // NOLINT
              "MUTK minor version must be less than 1000.");
static_assert(MUTK_VERSION_PATCH >= 0 && MUTK_VERSION_PATCH < 10000,  // NOLINT
              "MUTK patch version must be less than 10000.");

bool mutk::version_number_check_equal(int version_int) {
    return version_int == MUTK_VERSION_INTEGER;
}

TEST_CASE("version_number_check_equal") {
    CHECK(mutk::version_number_check_equal(MUTK_VERSION_INTEGER) == true);
    CHECK(mutk::version_number_check_equal(-1) == false);
}

int mutk::version_integer() { return MUTK_VERSION_INTEGER; }

TEST_CASE("version_integer") {
    CHECK(mutk::version_integer() == MUTK_VERSION_INTEGER);
}
