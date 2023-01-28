/*
# Copyright (c) 2023 Reed A. Cartwright <racartwright@gmail.com>
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

// NOTE: This header is in the lib directory because it is not part of the
// public API at this time.

#ifndef MUTK_JUNCTION_TREE_HPP
#define MUTK_JUNCTION_TREE_HPP

#include <mutk/graph.hpp>

namespace mutk::junction_tree {
struct component_t {
    std::vector<mutk::variable_t> variables;
    std::vector<float> edge_lengths;
};

using clique_t = std::vector<mutk::RelationshipGraph::vertex_descriptor>;

mutk::JunctionTree
create_junction_tree(const mutk::RelationshipGraph &graph,
    const std::vector<component_t> &components,
    const std::vector<clique_t> &elimination_order);

}

#endif // MUTK_JUNCTION_TREE_HPP
