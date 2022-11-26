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

#ifndef MUTK_INHERITANCE_MODEL_HPP
#define MUTK_INHERITANCE_MODEL_HPP

#include <type_traits>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace mutk {

class GraphBuilder;

class InheritanceModel {
 public:
    InheritanceModel();

    enum struct chromosome_type_t : int {};

    chromosome_type_t AddType(std::string_view name, int ploidy);

 private:
    std::vector<std::string> sexes_;
    std::unordered_map<std::string, chromosome_type_t> map_name_to_type_;

    std::vector<int> ploidies_;

    struct pattern_t {
        std::vector<chromosome_type_t> pattern;
        std::vector<int> discard;
    };

    std::vector<pattern_t> patterns_;

    friend GraphBuilder;
};

constexpr auto operator+(InheritanceModel::chromosome_type_t value) {
    return static_cast<std::underlying_type_t<InheritanceModel::chromosome_type_t>>(value);
}

} // namespace mutk

#endif