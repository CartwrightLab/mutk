/*
# Copyright (c) 2020 Reed A. Cartwright <reed@cartwright.ht>
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

/*

##PEDNG v0.1
##
#Indiv        Dad         Mom         Sex   Samples
A1@founder    .           .           1     .
A2            .           .           2     A2a A2b
A3            A1          A2          1     =
A4@gamete     A1:0.1      .           1     =
A5@clone      A3          .           1     =
A6@ploidy=1@founder .     .           1     =

*/

#ifndef MUTK_PEDIGREE_HPP
#define MUTK_PEDIGREE_HPP

#include <vector>
#include <string>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <optional>

#include <mutk/utility.hpp>

namespace mutk {

class Pedigree {
public:
    enum class Sex : char {
        Autosomal = 0, Male = 1, Female = 2,
        Invalid
    };

    static Sex parse_sex(const std::string &str);

    template<typename Range>
    static Pedigree parse_text(const Range &text);

    static Pedigree parse_table(const std::vector<std::vector<std::string>> &table);

    struct Member {
        std::string name;
        std::vector<std::string> tags;

        std::optional<std::string> dad;
        std::optional<double> dad_length;

        std::optional<std::string> mom;
        std::optional<double> mom_length;

        Sex sex;
 
        std::vector<std::string> samples;
    };

    using MemberTable = std::vector<Member>;

    Pedigree() = default;

    explicit Pedigree(MemberTable table) : table_{table} {
        for(MemberTable::size_type i=0; i < table_.size(); ++i) {
            auto ret = names_.emplace(table_[i].name, i);
            if(!ret.second) {
                throw std::invalid_argument("The name of a member of the pedigree is not unique: '" + table_[i].name + "'.");
            }
        }
    }

    std::size_t AddMember(Member member) {
        auto pos = table_.size();
        auto ret = names_.emplace(member.name, pos);
        if(!ret.second) {
            throw std::invalid_argument("The name of a member of the pedigree is not unique: '" + member.name + "'.");
        }
        table_.push_back(std::move(member));
        return pos;
    }

    std::size_t LookupMemberPosition(const std::string & child) const {
        auto it = names_.find(child);
        if(it == names_.end()) {
            return names_.size();
        }
        return it->second;
    }

    const Member* LookupMember(const std::string & child) const {
        auto it = names_.find(child);
        if(it == names_.end()) {
            return nullptr;
        }
        return &table_[it->second];
    }

    const Member& GetMember(std::size_t pos) const { return table_.at(pos); }

    std::size_t NumberOfMembers() const { return table_.size(); } 

    void Clear() {
        table_.clear();
        names_.clear();
    }

    const MemberTable& table() { return table_; }

private:
    MemberTable table_;
    std::unordered_map<std::string,MemberTable::size_type> names_;
};

inline
Pedigree::Sex Pedigree::parse_sex(const std::string &str) {
    static std::pair<std::string, Sex> keys[] = {
        {".", Sex::Invalid},
        {"0", Sex::Autosomal},
        {"1", Sex::Male},
        {"2", Sex::Female},
        {"male", Sex::Male},
        {"female", Sex::Female},
        {"autosomal", Sex::Autosomal}
    };
    return utility::key_switch_tuple(str, keys, keys[0]).second;
}

template<typename Range>
Pedigree Pedigree::parse_text(const Range &text) {
    using namespace boost;
    using namespace std;
    // Construct the tokenizer
    // token are separated by one or more <space>s or <tab>s
    // <newline>s end the row
    auto tokens = utility::make_tokenizer_dropempty(text, "\t ", "\n");

    auto token_it = tokens.begin();
    if(token_it == tokens.end() || *token_it != "##PEDNG") {
        throw std::invalid_argument("Pedigree parsing failed; "
            "unknown pedigree format; missing '##PEDNG' header line.");
    }

    // Work through tokens and build each row
    size_t k = 0;
    bool in_comment = false;
    vector<vector<string>> string_table;
    string_table.reserve(64);
    for(; token_it != tokens.end(); ++token_it) {
        const auto & token = *token_it;
        if(token == "\n") {
            k = 0;
            in_comment = false;
            continue;
        }
        if(in_comment) {
            // skip rows that are comments
            continue;
        }
        if(k == 0) {
            string_table.emplace_back();
            if(token[0] == '#') {
                in_comment = true;
                continue;
            }
            string_table.back().reserve(8);
        }
        // Add token to the current row
        string_table.back().push_back(token);
        k += 1;
    }

    return parse_table(string_table);
 }


} // namespace mutk

#endif //MUTK_PEDIGREE_HPP
