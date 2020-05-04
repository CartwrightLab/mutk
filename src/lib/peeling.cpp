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

void mutk::PeelUp::operator()(workspace_t *work, const TransitionMatrix &mat) const {
    assert(work != nullptr);
    work->lower[output_] = (mat * work->lower[a_].matrix()).array();
}

void mutk::PeelUpStar::operator()(workspace_t *work, const TransitionMatrix &mat) const {
    assert(work != nullptr);
    work->lower[output_] = (mat * work->lower[a_].matrix()).array() *
        work->lower[b_];
}

void mutk::PeelDown::operator()(workspace_t *work, const TransitionMatrix &mat) const {
    assert(work != nullptr);
    work->upper[output_] = (mat * (work->upper[a_] *
                         work->lower[b_]).matrix()).array();
}

void mutk::PeelStar::operator()(workspace_t *work, const TransitionMatrix &/*mat*/) const {
    assert(work != nullptr);
    work->upper[output_] = work->upper[a_] * work->lower[b_];
}

void mutk::PeelMate::operator()(workspace_t *work, const TransitionMatrix &/*mat*/) const {
    assert(work != nullptr);
    auto mate_width = work->upper[b_].size();
    auto pivot_width = work->lower[a_].size()/mate_width;

    Eigen::Map<TransitionMatrix> vmat{work->lower[a_].data(), pivot_width, mate_width};
    work->lower[output_] = (vmat * work->upper[b_].matrix()).array();
}
