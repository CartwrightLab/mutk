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

#include <string>
#include <filesystem>
#include <map>

#include <mutk/mutk.hpp>
#include <mutk/vcf.hpp>
#include <mutk/relationship_graph.hpp>

#include <CLI11.hpp>

#include "subcommand.hpp"

using namespace std::string_literals;

using InheritanceModel = mutk::InheritanceModel;

namespace {
struct args_t {
    double mu{1e-8};

    InheritanceModel chr_model{InheritanceModel::Autosomal};

    std::filesystem::path ped{};
    std::filesystem::path input{};
} args;
}  // anon namespace

int main(int argc, char *argv[]) {
    MUTK_RUNTIME_CHECK_VERSION_NUMBER_OR_RETURN();

    using namespace mutk::subcommand::string_literals;

    CLI::App app{"mutk modelfit v" MUTK_VERSION};

    #define ADD_OPTION_(name, desc) app.add_option(#name##_opt, args.name, desc, true)

    ADD_OPTION_(mu, "Germline mutation rate");
    ADD_OPTION_(ped, "Pedigree file");

    ADD_OPTION_(chr_model, "Chromosomal inheritance model")
        ->transform(CLI::CheckedTransformer(mutk::detail::CHR_MODEL_MAP, CLI::ignore_case));

    app.add_option("input", args.input, "Input file");
    #undef ADD_OPTION_

    CLI11_PARSE(app, argc, argv);

    mutk::BcfReader reader{args.input};

    auto samples = reader.samples();

    std::vector<const char*> known_samples(samples.first, samples.first+samples.second);

    auto pedigree = mutk::Pedigree::parse_file(args.ped);

    mutk::RelationshipGraph graph;
    graph.ConstructGraph(pedigree, known_samples, args.chr_model, args.mu, args.mu, false);

    graph.PrintGraph(std::cout);

    return EXIT_SUCCESS;
}
