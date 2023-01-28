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

#ifndef MUTK_GRAPH_BUILDER_HPP
#define MUTK_GRAPH_BUILDER_HPP

#include <string>
#include <vector>
#include <unordered_map>

#include <boost/graph/adjacency_list.hpp>

#include "inheritance_model.hpp"
#include "graph.hpp"
#include "potential.hpp"


namespace mutk {

class GraphBuilder {
 public:
    GraphBuilder();

    // Add families as biconnected components
    member_id_t AddSingle(const std::string &name, const std::string &sex, const std::vector<std::string> &sample_names);

    member_id_t AddPair(const std::string &name, const std::string &sex, const std::vector<std::string> &sample_names,
        const std::string &parent_name, float mutation_scale);

    member_id_t AddTrio(const std::string &name, const std::string &sex, const std::vector<std::string> &sample_names,
        const std::string &parent_name_a, float mutation_scale_a,
        const std::string &parent_name_b, float mutation_scale_b);

    // Set known samples
    void SetSamples(const std::vector<std::string> &sample_names);

    // Build the final relationship graph based on an inhertiance model
    RelationshipGraph BuildGraph(const InheritanceModel &model, float mu);

 private:
    member_id_t LookupName(const std::string &name);

    RelationshipGraph CreateInitialGraph(const InheritanceModel &model, float mu);

    struct family_t {
        std::vector<member_id_t> members;
        std::vector<float> scales;
    };

    int member_counter_{0};

    std::unordered_map<std::string, member_id_t> map_member_name_to_id_;
    std::vector<std::string> member_names_;
    std::vector<std::string> member_sexes_;
    std::vector<std::vector<std::string>> member_input_samples_;
    std::vector<std::vector<sample_id_t>> member_data_samples_;

    std::unordered_map<std::string, sample_id_t> map_sample_name_to_id_;
    std::vector<std::string> sample_names_;

    std::vector<family_t> families_;
};

}

#endif // MUTK_GRAPH_BUILDER_HPP
