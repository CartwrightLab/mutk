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

#include <mutk/mutk.hpp>
#include <mutk/vcf.hpp>
#include <mutk/relationship_graph.hpp>
#include <mutk/memory.hpp>
#include <mutk/utility.hpp>

#include <CLI11.hpp>

#include "subcommand.hpp"

using namespace std::string_literals;

using InheritanceModel = mutk::InheritanceModel;
using Potential = mutk::RelationshipGraph::Potential;

namespace {
struct args_t {
    double mu{1e-8};

    double theta{0.001};
    double ref_bias_hom{0.0};
    double ref_bias_het{0.0};
    double ref_bias_hap{0.0};
    
    double seq_error{0.005};
    double seq_bias{0.0};
    double seq_overdisp_hom{0.0};
    double seq_overdisp_het{0.0};

    InheritanceModel chr_model{InheritanceModel::Autosomal};

    std::filesystem::path ped{};
    std::filesystem::path output{};
    std::filesystem::path input{};
};
}  // anon namespace

int main(int argc, char *argv[]) {
    MUTK_RUNTIME_CHECK_VERSION_NUMBER_OR_RETURN();

    using namespace mutk::subcommand::string_literals;

    args_t args;

    CLI::App app{"mutk modelfit v" MUTK_VERSION};

    #define ADD_OPTION_(name, desc) app.add_option(#name##_opt, args.name, desc, true)

    ADD_OPTION_(mu, "Germline mutation rate");

    ADD_OPTION_(theta, "Population diversity");

    ADD_OPTION_(ref_bias_hom, "Ascertainment bias for reference homozygotes");
    ADD_OPTION_(ref_bias_het, "Ascertainment bias for reference heterozygotes");
    ADD_OPTION_(ref_bias_hap, "Ascertainment bias for reference haploids");

    // ADD_OPTION_(seq_error, "Sequencing error rate (per base-call)");
    // ADD_OPTION_(seq_bias, "Sequencing reference bias");
    // ADD_OPTION_(seq_overdisp_hom, "Sequencing overdispersion for homozygotes");
    // ADD_OPTION_(seq_overdisp_het, "Sequencing overdispersion for heterozygotes");

    ADD_OPTION_(ped, "Pedigree file");
    ADD_OPTION_(chr_model, "Chromosomal inheritance model")
        ->transform(CLI::CheckedTransformer(mutk::detail::CHR_MODEL_MAP, CLI::ignore_case));

    ADD_OPTION_(output, "Output file");

    #undef ADD_OPTION_

    app.add_option("input", args.input, "Input file");

    CLI11_PARSE(app, argc, argv);

    mutk::vcf::Reader reader{args.input};

    auto samples = reader.samples();

    std::vector<const char*> known_samples(samples.first, samples.first+samples.second);

    auto pedigree = mutk::Pedigree::parse_file(args.ped);

    mutk::RelationshipGraph graph;
    graph.ConstructGraph(pedigree, known_samples, args.chr_model, args.mu, args.mu, false);

    graph.ConstructPeeler();

    int iret = reader.SetSamples(graph.SampleNames());
    assert(iret == 0);

    // allocate space for PL data
    int num_samples = reader.samples().second;
    int n_pl_capacity = 15*num_samples;
    auto pl_buf = mutk::vcf::make_buffer<int>(n_pl_capacity);

    // allocate workspace
    auto work = graph.CreateWorkspace();

    mutk::mutation::KAllelesModel model(5.0, args.theta,
        args.ref_bias_hom, args.ref_bias_het, args.ref_bias_hap);

    reader([&](const bcf_hdr_t *header, bcf1_t *record) {
        if(record->n_allele > 5) {
            // we currently do not support locations with more than
            // five alleles
            return;
        }
        bcf_unpack(record, BCF_UN_ALL);

        int n_pl = mutk::vcf::get_format_int32(header, record, "PL", &pl_buf);
        if(n_pl <= 0) {
            // PL tag is missing, so we do nothing at this time
            return;
        }
        const int num_genotypes = n_pl / num_samples;
        assert(n_pl % num_samples == 0);
        
        auto pl = mutk::wrap_tensor(pl_buf.data.get(), num_genotypes, num_samples);

        mutk::tensor_index_t haploid_sz = mutk::dim_width<1>(record->n_allele);
        mutk::tensor_index_t diploid_sz = mutk::dim_width<2>(record->n_allele);

        auto founder1 = model.CreatePriorHaploid(haploid_sz);
        auto founder2 = model.CreatePriorDiploid(haploid_sz);

        work.widths = {1, haploid_sz, diploid_sz};

        using mutk::utility::unphredf;
        for(int i=0; i < graph.potentials().size(); ++i) {
            const auto &pot = graph.potentials()[i];
            if(pot.type == Potential::LikelihoodDiploid) {
                work.stack[i].resize(diploid_sz);
                // check for missing data
                if(mutk::vcf::is_missing(pl(i,0))) {
                    // If PLs are missing for this site, set everything to 1.
                    work.stack[i].setConstant(1.0f);
                } else {                
                    int width = 0;
                    for(;width < num_genotypes; ++width) {
                        if(mutk::vcf::is_vector_end(pl(i,width))) {
                            break;
                        }
                    }
                    if(width != diploid_sz) {
                        // PL tag is not wide enough, we will skip the site
                        return;
                    }
                    for(int k=0; k < diploid_sz; ++k) {
                        assert(!mutk::vcf::is_missing(pl(i,k)));
                        work.stack[i](k) = unphredf(pl(i,k));
                    }
                }
            } else if(pot.type == Potential::LikelihoodHaploid) {
                work.stack[i].resize(haploid_sz);
                // check for missing data
                if(mutk::vcf::is_missing(pl(i,0))) {
                    // If PLs are missing for this site, set everything to 1.
                    work.stack[i].setConstant(1.0f);
                } else {
                    // Measure the number of genotype of this
                    int width = 0;
                    for(;width < num_genotypes; ++width) {
                        if(mutk::vcf::is_vector_end(pl(i,width))) {
                            break;
                        }
                    }
                    if(width == haploid_sz) {
                        // haploid PLs are encoded directly
                        for(int k=0; k < haploid_sz; ++k) {
                            work.stack[i](k) = unphredf(pl(i,k));
                        }
                    } else if(width == diploid_sz) {
                        // haploid PLs are encoded as homozygous diploids
                        int m = 0;
                        for(int k=0; k < haploid_sz; ++k) {
                            work.stack[i](k) = unphredf(pl(i,m));
                            m += k+2;
                        }
                    } else {
                        // PL tag is not wide enough, we will skip the site
                        return;
                    }
                }
            } else if(pot.type == Potential::FounderDiploid) {
                work.stack[i] = founder2;
            } else if(pot.type == Potential::FounderHaploid) {
                work.stack[i] = founder1;
            } else {
                work.stack[i] = model.CreatePotential(haploid_sz, pot, mutk::mutation::ANY);
            }
        }
        float value = graph.PeelForward(&work);
        std::cout << record->rid << " " << record->pos << " " << value << "\n";
    });

    // Go thorough the likelihood potentials and fill them with data from PL
    //    (1) If the likelihood is haploid, will need to check if it is encoded as a diploid.
    //    (2) If first element of a column is missing, fill the liklehood with 1s.
    //    (3) Convert PLs to normalized probabilities.

    // Initialize the priors based on n.
    // Initialize the mutation matrices based on n.

    // Store dimensions of unit matrices
    // Store founders and mutation matrices by value of n.
    // reshape potentials by elimination order

    {
        auto p = model.CreatePriorDiploid(2);
        mutk::Tensor<2> m1 = model.CreateMatrix(2, 1e-8f, mutk::mutation::ANY);
        mutk::Tensor<2> m2 = model.CreateMatrix(2, 2e-8f, mutk::mutation::ANY);
        mutk::detail::potential_t pot(Potential::ChildDiploidDiploid, 0, 1, 2e-8f, 2, 2e-8f);
        pot.shuffle = {0,1,2};
        mutk::Tensor<3> m3 = model.CreatePotential(2, pot, mutk::mutation::ANY).
            reshape(mutk::tensor_dims(3,3,3));
        double val = 0.0;
        for(int d = 0; d < 3; ++d) {
            for(int m = 0; m < 3; ++m) {
                double d1 = 0.0;
                for(int c = 0; c < 3; ++c) {
                    d1 += m3(c,d,m);
                }
                double d2 = 0.0;
                for(int c = 0; c < 3; ++c) {
                    d2 += m1(c,d);
                }
                double d3 = 0.0;
                for(int c = 0; c < 3; ++c) {
                    d3 += m1(c,m);
                }
                val += d1*d2*d3*p(d)*p(m);
            }
        }
        std::cout << val << " " << log(val) << std::endl;
    }

    return EXIT_SUCCESS;
}
