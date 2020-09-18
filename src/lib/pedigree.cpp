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
#include <doctest/doctest.h>

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
            // LCOV_EXCL_START
            if(tokens.begin() == tokens.end()) {
                throw std::invalid_argument("Pedigree parsing failed. Row "
                    + std::to_string(row_num) + " has empty child name."
                );
            }
            // LCOV_EXCL_STOP
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
        // We defer percent decoding to the newick parser
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

// LCOV_EXCL_START
TEST_CASE("[libmutk] Pedigree::parse_sex") {
    CHECK(Pedigree::parse_sex(".") == Pedigree::Sex::Invalid);
    CHECK(Pedigree::parse_sex("0") == Pedigree::Sex::Autosomal);
    CHECK(Pedigree::parse_sex("1") == Pedigree::Sex::Male);
    CHECK(Pedigree::parse_sex("2") == Pedigree::Sex::Female);
    CHECK(Pedigree::parse_sex("autosomal") == Pedigree::Sex::Autosomal);
    CHECK(Pedigree::parse_sex("male") == Pedigree::Sex::Male);
    CHECK(Pedigree::parse_sex("female") == Pedigree::Sex::Female);
    CHECK(Pedigree::parse_sex("a") == Pedigree::Sex::Autosomal);
    CHECK(Pedigree::parse_sex("m") == Pedigree::Sex::Male);
    CHECK(Pedigree::parse_sex("f") == Pedigree::Sex::Female);

    CHECK(Pedigree::parse_sex("malewrong") == Pedigree::Sex::Invalid);
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
TEST_CASE("[libmutk] Pedigree::parse_text") {
    const char ped[] = 
        "##PEDNG v1.0\n"
        "A . . 1 .\n"
        "B@founder\t.:0.1\t.\t2\tB1\n"
        "C    A    B    1    C1    C2\n"
        "D A:0.01 B:0.5 2\t=\n"
        "E@founder@haploid . . 2 =\n"
        "%46 %41 %42 1 %46\n"
        "%40@diploid %3a:0.1 %2e 1 %3D\n"
    ;
    Pedigree pedigree;
    REQUIRE_NOTHROW(pedigree = Pedigree::parse_text(ped));
    REQUIRE(pedigree.NumberOfMembers() == 7);
    CHECK(pedigree.table().data() == &pedigree.GetMember(0));

    CHECK(pedigree.GetMember(0).name == "A");
    CHECK(pedigree.GetMember(0).tags.empty());
    CHECK(pedigree.GetMember(0).dad == std::nullopt);
    CHECK(pedigree.GetMember(0).dad_length == std::nullopt);
    CHECK(pedigree.GetMember(0).mom == std::nullopt);
    CHECK(pedigree.GetMember(0).mom_length == std::nullopt);
    CHECK(pedigree.GetMember(0).sex == Pedigree::Sex::Male);
    CHECK(pedigree.GetMember(0).samples.empty());

    CHECK(pedigree.GetMember(1).name == "B");
    REQUIRE(pedigree.GetMember(1).tags.size() == 1);
    CHECK(pedigree.GetMember(1).tags[0] == "founder");
    CHECK(pedigree.GetMember(1).dad == std::nullopt);
    CHECK(pedigree.GetMember(1).dad_length == std::nullopt);
    CHECK(pedigree.GetMember(1).mom == std::nullopt);
    CHECK(pedigree.GetMember(1).mom_length == std::nullopt);
    CHECK(pedigree.GetMember(1).sex == Pedigree::Sex::Female);
    REQUIRE(pedigree.GetMember(1).samples.size() == 1);
    CHECK(pedigree.GetMember(1).samples[0] == "B1");

    CHECK(pedigree.GetMember(2).name == "C");
    CHECK(pedigree.GetMember(2).tags.empty());
    CHECK(pedigree.GetMember(2).dad == "A");
    CHECK(pedigree.GetMember(2).dad_length == std::nullopt);
    CHECK(pedigree.GetMember(2).mom == "B");
    CHECK(pedigree.GetMember(2).mom_length == std::nullopt);
    CHECK(pedigree.GetMember(2).sex == Pedigree::Sex::Male);
    REQUIRE(pedigree.GetMember(2).samples.size() == 2);
    CHECK(pedigree.GetMember(2).samples[0] == "C1");
    CHECK(pedigree.GetMember(2).samples[1] == "C2");

    CHECK(pedigree.GetMember(3).name == "D");
    CHECK(pedigree.GetMember(3).tags.empty());
    CHECK(pedigree.GetMember(3).dad == "A");
    CHECK(pedigree.GetMember(3).dad_length == 0.01);
    CHECK(pedigree.GetMember(3).mom == "B");
    CHECK(pedigree.GetMember(3).mom_length == 0.5);
    CHECK(pedigree.GetMember(3).sex == Pedigree::Sex::Female);
    REQUIRE(pedigree.GetMember(3).samples.size() == 1);
    CHECK(pedigree.GetMember(3).samples[0] == "D");

    CHECK(pedigree.GetMember(4).name == "E");
    REQUIRE(pedigree.GetMember(4).tags.size() == 2);
    CHECK(pedigree.GetMember(4).tags[0] == "founder");
    CHECK(pedigree.GetMember(4).tags[1] == "haploid");
    CHECK(pedigree.GetMember(4).dad == std::nullopt);
    CHECK(pedigree.GetMember(4).dad_length == std::nullopt);
    CHECK(pedigree.GetMember(4).mom == std::nullopt);
    CHECK(pedigree.GetMember(4).mom_length == std::nullopt);
    CHECK(pedigree.GetMember(4).sex == Pedigree::Sex::Female);
    REQUIRE(pedigree.GetMember(4).samples.size() == 1);
    CHECK(pedigree.GetMember(4).samples[0] == "E");

    CHECK(pedigree.GetMember(5).name == "F");
    CHECK(pedigree.GetMember(5).tags.empty());
    CHECK(pedigree.GetMember(5).dad == "A");
    CHECK(pedigree.GetMember(5).dad_length == std::nullopt);
    CHECK(pedigree.GetMember(5).mom == "B");
    CHECK(pedigree.GetMember(5).mom_length == std::nullopt);
    CHECK(pedigree.GetMember(5).sex == Pedigree::Sex::Male);
    REQUIRE(pedigree.GetMember(5).samples.size() == 1);
    CHECK(pedigree.GetMember(5).samples[0] == "%46");    

    CHECK(pedigree.GetMember(6).name == "@");
    REQUIRE(pedigree.GetMember(6).tags.size() == 1);
    CHECK(pedigree.GetMember(6).tags[0] == "diploid");
    CHECK(pedigree.GetMember(6).dad == ":");
    CHECK(pedigree.GetMember(6).dad_length == 0.1);
    CHECK(pedigree.GetMember(6).mom == ".");
    CHECK(pedigree.GetMember(6).mom_length == std::nullopt);
    CHECK(pedigree.GetMember(6).sex == Pedigree::Sex::Male);
    REQUIRE(pedigree.GetMember(6).samples.size() == 1);
    CHECK(pedigree.GetMember(6).samples[0] == "%3D");

    CHECK_THROWS_AS(Pedigree::parse_text(""), std::invalid_argument);
    CHECK_THROWS_AS(Pedigree::parse_text("#PEDNG"), std::invalid_argument);
    CHECK_THROWS_AS(Pedigree::parse_text("##PEDNG v1.0\nA\t.\t."), std::invalid_argument);
    CHECK_THROWS_AS(Pedigree::parse_text("##PEDNG v1.0\n@\t.\t.\t1\t."), std::invalid_argument);
    CHECK_THROWS_AS(Pedigree::parse_text("##PEDNG v1.0\nA\t:0.1\t.\t1\t."), std::invalid_argument);
    CHECK_THROWS_AS(Pedigree::parse_text("##PEDNG v1.0\nA\tB:q\t.\t1\t."), std::invalid_argument);
    CHECK_THROWS_AS(Pedigree::parse_text("##PEDNG v1.0\nA\t.\t:0.1\t1\t."), std::invalid_argument);
    CHECK_THROWS_AS(Pedigree::parse_text("##PEDNG v1.0\nA\t.\tB:q\t1\t."), std::invalid_argument);
    CHECK_THROWS_AS(Pedigree::parse_text("##PEDNG v1.0\nA\t.\t.\t.\t."), std::invalid_argument);

}
// LCOV_EXCL_STOP

} // namespace mutk
