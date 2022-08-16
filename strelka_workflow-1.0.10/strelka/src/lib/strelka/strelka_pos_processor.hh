// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///

///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#ifndef __STRELKA_POS_PROCESSOR_HH
#define __STRELKA_POS_PROCESSOR_HH

#include "strelka_sample_type.hh"
#include "strelka_streams.hh"

#include "starling_common/starling_pos_processor_base.hh"
#include "strelka/strelka_shared.hh"


///
///
struct strelka_pos_processor : public starling_pos_processor_base {

    typedef starling_pos_processor_base base_t;

    strelka_pos_processor(const strelka_options& opt,
                          const strelka_deriv_options& dopt,
                          const reference_contig_segment& ref,
                          const strelka_streams& client_io);

    ~strelka_pos_processor();

private:
    void
    process_pos_variants(const pos_t pos) {
        process_pos_indel_single_sample(pos,STRELKA_SAMPLE_TYPE::NORMAL);
        process_pos_indel_single_sample(pos,STRELKA_SAMPLE_TYPE::TUMOR);

        process_pos_snp_single_sample(pos,STRELKA_SAMPLE_TYPE::NORMAL);
        process_pos_snp_single_sample(pos,STRELKA_SAMPLE_TYPE::TUMOR);

        process_pos_snp_somatic(pos);
        process_pos_indel_somatic(pos);
    }

    void
    process_pos_snp_somatic(const pos_t pos);

    void
    process_pos_indel_somatic(const pos_t pos);

    void
    write_counts(const pos_range& output_report_range) const;

    /////////////////////////////

    // keep some of the original strelka classes handy so we don't
    // have to down-cast:
    const strelka_options& _opt;
    const strelka_deriv_options& _dopt;
    const strelka_streams& _client_io;

    extra_position_data _tier2_epd[MAX_SAMPLE];
};


#endif
