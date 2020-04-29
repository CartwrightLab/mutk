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

#include <mutk/pedigree.hpp>

namespace mutk {

Pedigree Pedigree::parse_table(const std::vector<std::vector<std::string>> &table) {
    using utility::percent_decode;

   // Build Pedigree Object
    Pedigree ret;
    int row_num = 0;
    for(auto &&row : table) {
        Member member;
        row_num += 1;
        if(row.empty()) {
            continue;
        }
        if(row.size() < 5) {
            throw std::invalid_argument("Pedigree parsing failed. Row "
                + std::to_string(row_num) + " has "
                + std::to_string(row.size()) + " column(s) instead of 5 or more columns."
            );
        }
        // separate tags from member name
        {
            if(row[0][0] == '@') {
                throw std::invalid_argument("Pedigree parsing failed. Row "
                    + std::to_string(row_num) + " has empty child name."
                );                
            }
            auto tokens = utility::make_tokenizer_dropempty(row[0], "@", "");
            if(tokens.begin() == tokens.end()) {
                throw std::invalid_argument("Pedigree parsing failed. Row "
                    + std::to_string(row_num) + " has empty child name."
                );
            }
            auto it = tokens.begin();
            member.name = percent_decode(*(it++));
            for(;it != tokens.end();++it) {
                member.tags.push_back(percent_decode(*it));
            }
        }
        // process dad column
        {
            auto pos = row[1].find(':');
            member.dad = row[1].substr(0,pos);
            if(member.dad->empty()) {
                throw std::invalid_argument("Pedigree parsing failed. Row "
                    + std::to_string(row_num) + " has empty dad name."
                );
            }
            if(*member.dad == ".") {
                member.dad = std::nullopt;
            } else {
                member.dad = percent_decode(*member.dad);
                if(pos != std::string::npos) {
                    char *str_end;
                    member.dad_length = std::strtod(row[1].c_str()+pos+1, &str_end);
                    if(pos+1 == row[1].size() || str_end != row[1].c_str()+row[1].size()) {
                        throw std::invalid_argument("Pedigree parsing failed. Row "
                            + std::to_string(row_num) + " has invalid dad length."
                        );
                    }
                }
            }
        }
        // process mom column
        {
            auto pos = row[2].find(':');
            member.mom = row[2].substr(0,pos);
            if(member.mom->empty()) {
                throw std::invalid_argument("Pedigree parsing failed. Row "
                    + std::to_string(row_num) + " has empty mom name."
                );
            }
            if(*member.mom == ".") {
                member.mom = std::nullopt;
            } else {
                member.mom = percent_decode(*member.mom);
                if(pos != std::string::npos) {
                    char *str_end;
                    member.mom_length = std::strtod(row[2].c_str()+pos+1, &str_end);
                    if(pos+1 == row[2].size() || str_end != row[2].c_str()+row[2].size()) {
                        throw std::invalid_argument("Pedigree parsing failed. Row "
                            + std::to_string(row_num) + " has invalid mom length."
                        );
                    }
                }
            }
        }
        // process chromosmal sex column
        {
            member.sex = parse_sex(row[3]);
            if(member.sex == Sex::Invalid) {
                throw std::invalid_argument("Pedigree parsing failed. Row "
                    + std::to_string(row_num) + " has invalid sex."
                );
            }
        }
        // Process samples
        for(auto it = row.begin()+4; it != row.end();++it) {
             if(*it == "=") {
                member.samples.push_back(member.name);
            } else if(*it != ".") {
                member.samples.push_back(*it);
            }
        }
        ret.AddMember(member);
    }
    return ret;
}

} // namespace mutk
