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

#ifndef __GVCF_LOCUS_INFO_HH
#define __GVCF_LOCUS_INFO_HH


#include "blt_common/position_snp_call_pprob_digt.hh"
#include "starling_common/align_path.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "starling_common/starling_shared.hh"

#include <bitset>
#include <iosfwd>


namespace VCF_FILTERS {

enum index_t {
    IndelConflict,
    SiteConflict,
    LowGQX,
    HighBaseFilt,
    HighDepth,
    HighSNVSB,
    HighSNVHPOL,
    HighRefRep,
    SIZE
};

inline
const char*
get_label(const unsigned idx) {
    switch (idx) {
    case HighDepth: return "HighDepth";
    case LowGQX: return "LowGQX";
    case HighSNVSB: return "HighSNVSB";
    case HighSNVHPOL: return "HighSNVHPOL";
    case HighBaseFilt: return "HighDPFRatio";
    case HighRefRep: return "HighREFREP";
    case IndelConflict: return "IndelConflict";
    case SiteConflict: return "SiteConflict";
    default:
        assert(0);
        return NULL;
    }
}
}



struct shared_modifiers {

    shared_modifiers() { clear(); }

    void
    set_filter(const VCF_FILTERS::index_t i) {
        filters.set(i);
    }

    void
    write_filters(std::ostream& os) const;

    void
    clear() {
        filters.reset();
    }

    int gqx;
    int gq;
    unsigned max_gt;

    std::bitset<VCF_FILTERS::SIZE> filters;
};


std::ostream& operator<<(std::ostream& os,const shared_modifiers& shmod);


struct indel_modifiers : public shared_modifiers {
    indel_modifiers() { clear(); }

    void
    clear() {
        shared_modifiers::clear();
        is_overlap=false;
        ploidy.clear();
    }

    ALIGNPATH::path_t cigar;

    bool is_overlap;
    std::vector<unsigned> ploidy;
};



namespace MODIFIED_SITE_GT {

enum index_t {
    NONE,
    UNKNOWN,
    ZERO,
    ONE
};

inline
const char*
get_label(const unsigned idx) {
    switch (static_cast<index_t>(idx)) {
    case ZERO: return "0";
    case ONE: return "1";
    case UNKNOWN: return ".";
    default:
        assert(0);
        return NULL;
    }
}
}

struct site_modifiers : public shared_modifiers {

    site_modifiers() { clear(); }

    void
    clear() {
        shared_modifiers::clear();
        is_unknown=true;
        is_covered=false;
        is_used_covered=false;
        is_zero_ploidy=false;
        is_block=false;
        modified_gt=MODIFIED_SITE_GT::NONE;
    }

    bool
    is_gqx() const {
        return ((!is_unknown) && is_used_covered && (!is_zero_ploidy));
    }

    bool is_unknown;
    bool is_covered;
    bool is_used_covered;
    bool is_zero_ploidy;
    bool is_block;

    MODIFIED_SITE_GT::index_t modified_gt;
};


std::ostream& operator<<(std::ostream& os,const site_modifiers& smod);



struct indel_info {

    void
    init(const pos_t init_pos,
         const indel_key& init_ik,
         const starling_diploid_indel_core& init_dindel,
         const starling_indel_report_info& init_iri,
         const starling_indel_sample_report_info& init_isri)
    {
        pos=(init_pos);
        ik=(init_ik);
        dindel=(init_dindel);
        iri=(init_iri);
        isri=(init_isri);
//        MQ=0.0;
//        ReadPosRankSum=0.0;
//        BaseQRankSum=0.0;
//        MQRankSum=0.0;
        imod.clear();
        MQ = 0.0;
        ReadPosRankSum = 0.0;
        BaseQRankSum = 0.0;
        MQRankSum = 0.0;
    }

    const char*
    get_gt() {
        if (imod.is_overlap) {
            return "1/2";
        }
        return STAR_DIINDEL::get_gt_label(imod.max_gt);
    }

    // the site ploidy within the indel at offset x
    unsigned
    get_ploidy(const unsigned offset) {
        if (! imod.is_overlap) {
            using namespace STAR_DIINDEL;
            switch (dindel.max_gt) {
            case HOM: return 0;
            case HET: return 1;
            case NOINDEL: return 2;
            }
            assert(0);
        } else {
            assert(offset<imod.ploidy.size());
            return imod.ploidy[offset];
        }
        return 2;
    }


    pos_t pos;
    indel_key ik;
    starling_diploid_indel_core dindel;
    starling_indel_report_info iri;
    starling_indel_sample_report_info isri;
//    double MQ;               // RMS of mapping qualities
//    double ReadPosRankSum;   // Uses Mann-Whitney Rank Sum Test for the distance from the end of the read containing an alternate allele.
//    double BaseQRankSum;     // Uses Mann-Whitney Rank Sum Test for BQs (ref bases vs alternate alleles)
//    double MQRankSum;        // Uses Mann-Whitney Rank Sum Test for MQs (ref bases vs alternate alleles)
    indel_modifiers imod;
    double MQ;               // RMS of mapping qualities

    //only meaningful for het calls
    double ReadPosRankSum;   // Uses Mann-Whitney Rank Sum Test for the distance from the end of the read containing an alternate allele.
    double BaseQRankSum;     // Uses Mann-Whitney Rank Sum Test for BQs (ref bases vs alternate alleles)
    double MQRankSum;

};


//Data structure defining parameters for a single site to be used for writing in gvcf_aggregator
struct site_info {

    site_info()
        : pos(0)
        , ref('N')
        , n_used_calls(0)
        , n_unused_calls(0)
        , hpol(0)
        , hapscore(0)
        , MQ(0)
        , ReadPosRankSum(0)
        , BaseQRankSum(0)
        , MQRankSum(0)
    {
        for (unsigned i(0); i<N_BASE; ++i) known_counts[i] = 0;
    }

    void
    init(const pos_t init_pos,
         const char init_ref,
         const snp_pos_info& good_pi,
         const bool used_allele_count_min_qscore) {

        pos=(init_pos);
        ref=(init_ref);
        good_pi.get_known_counts(known_counts,used_allele_count_min_qscore);
        smod.clear();
    }


    const char*
    get_gt() const {
        if       (smod.modified_gt != MODIFIED_SITE_GT::NONE) {
            return MODIFIED_SITE_GT::get_label(smod.modified_gt);
        } else if (smod.is_unknown || (!smod.is_used_covered)) {
            return ".";
        } else {
            unsigned print_gt(smod.max_gt);
            if (smod.is_block) {
                print_gt = dgt.ref_gt;
            }
            return DIGT::get_vcf_gt(print_gt,dgt.ref_gt);
        }
    }

    bool
    is_qual() const {
        return ((!smod.is_block) && (!smod.is_unknown) && smod.is_used_covered && (!smod.is_zero_ploidy) && (dgt.ref_gt != smod.max_gt));
    }


    pos_t pos;
    char ref;
    unsigned n_used_calls;
    unsigned n_unused_calls;
    boost::array<unsigned,N_BASE> known_counts;
    diploid_genotype dgt;
    unsigned hpol;
    double hapscore;
    double MQ;				 // RMS of mapping qualities

    //only meaningful for het calls
    double ReadPosRankSum;   // Uses Mann-Whitney Rank Sum Test for the distance from the end of the read containing an alternate allele.
    double BaseQRankSum;     // Uses Mann-Whitney Rank Sum Test for BQs (ref bases vs alternate alleles)
    double MQRankSum;        // Uses Mann-Whitney Rank Sum Test for MQs (ref bases vs alternate alleles)

    site_modifiers smod;
};

#endif
