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

#include <mutk/peeling.hpp>

using workspace_t = mutk::workspace_t;

void mutk::PeelUp::operator()(workspace_t *work) const {
    assert(work != nullptr);
    work->lower[output_] = (mat_ * work->lower[a_].matrix()).array();
}

void mutk::PeelUpStar::operator()(workspace_t *work) const {
    assert(work != nullptr);
    work->lower[output_] = (mat_ * work->lower[a_].matrix()).array() *
        work->lower[b_];
}


void mutk::PeelDown::operator()(workspace_t *work) const {
    assert(work != nullptr);
    work->upper[output_] = (mat_ * (work->upper[a_] *
                         work->lower[b_]).matrix()).array();
}

void mutk::PeelLeft::Operator()(workspace_t *work) const {
    assert(work != nullptr);
    

}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_father(workspace_t &work, const family_members_t &family,
                          const TransitionMatrixVector &mat) {
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.temp_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    // Include Mom
    auto mom_width = work.upper[mom].size();
    auto dad_width = work.temp_buffer.size()/mom_width;
    work.temp_buffer.resize(mom_width, dad_width);
    work.lower[dad] *= (work.temp_buffer.matrix().transpose() * (work.upper[mom] *
                        work.lower[mom]).matrix()).array();

    //work.temp_buffer.resize(work.temp_buffer.size(), 1);
}


// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_mother(workspace_t &work, const family_members_t &family,
                          const TransitionMatrixVector &mat) {
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.temp_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(size_t i = 3; i < family.size(); ++i) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    // Include Mom
    auto dad_width = work.upper[dad].size();
    auto mom_width = work.temp_buffer.size()/dad_width;
    work.temp_buffer.resize(mom_width, dad_width);
    work.lower[mom] *= (work.temp_buffer.matrix() * (work.upper[dad] *
                        work.lower[dad]).matrix()).array();
    //work.temp_buffer.resize(work.temp_buffer.size(), 1);
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_mother_fast(workspace_t &work,
                               const family_members_t &family, const TransitionMatrixVector &mat) {
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.temp_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    // Include Dad
    auto dad_width = work.upper[dad].size();
    auto mom_width = work.temp_buffer.size()/dad_width;
    work.temp_buffer.resize(mom_width, dad_width);
    work.lower[mom] = (work.temp_buffer.matrix() * (work.upper[dad] *
                       work.lower[dad]).matrix()).array();
    //work.temp_buffer.resize(work.temp_buffer.size(), 1);
}

// Family Order: Father, Mother, Child, Child2, ....
void dng::peel::to_child(workspace_t &work, const family_members_t &family,
                         const TransitionMatrixVector &mat) {
    assert(family.size() >= 4);
    auto dad = family[0];
    auto mom = family[1];
    auto child = family[2];
    // Parents
    work.temp_buffer = kroneckerProduct((work.lower[dad] *
                                           work.upper[dad]).matrix(),
                                          (work.lower[mom] * work.upper[mom]).matrix()).array();
    // Sum over fullsibs
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.temp_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }

    work.upper[child] = (mat[child].transpose() *
                         work.temp_buffer.matrix()).array();
}

// Family Order: Father, Mother, CHild
void dng::peel::to_child_fast(workspace_t &work, const family_members_t &family,
                              const TransitionMatrixVector &mat) {
    assert(family.size() == 3);
    auto dad = family[0];
    auto mom = family[1];
    auto child = family[2];
    work.upper[child] = (mat[child].transpose() * kroneckerProduct(
                             (work.lower[dad] * work.upper[dad]).matrix(),
                             (work.lower[mom] * work.upper[mom]).matrix())).array();
}


