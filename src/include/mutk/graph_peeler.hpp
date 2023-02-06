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

#ifndef MUTK_RELATIONSHIP_GRAPH_HPP
#define MUTK_RELATIONSHIP_GRAPH_HPP

#include "message.hpp"
#include "graph.hpp"

#include <cmath>
#include <string>
#include <vector>

namespace mutk {

struct workspace_t {
    std::vector<mutk::message_t> messages;
};

class GraphPeeler {
public:
    GraphPeeler() = default;

    static GraphPeeler Create(RelationshipGraph graph);

    float PeelForward(workspace_t &work) const;

    template<class Arg>
    void SetModelPotentials(workspace_t &work, message_size_t n, Arg arg) const;

    void SetDataPotentials(workspace_t &work, message_size_t n,
        const std::vector<mutk::message_t> &data) const;

    workspace_t CreateWorkspace() const;

    // const std::vector<potential_t>& potentials() const {
    //     return potentials_;
    // }

    const auto & graph() const {
        return graph_;
    }

    const auto & junction_tree() const {
        return tree_;
    }

protected:
    RelationshipGraph graph_;
    JunctionTree tree_;

private:

};

} // namespace mutk

#endif // MUTK_RELATIONSHIP_GRAPH_HPP
