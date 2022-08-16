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
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#include "blt_common/blt_arg_parse_util.hh"
#include "blt_common/blt_arg_validate.hh"
#include "starling_common/starling_shared.hh"


void
legacy_starling_arg_parse(arg_data& ad,
                          starling_options& opt) {

    const prog_info& pinfo(ad.pinfo);

    if ((ad.size()==1) ||
        ((ad.size()==2) && ((ad.argstr[1] == "-h") || (ad.argstr[1] == "-help") ||
                            (ad.argstr[1] == "--help") || (ad.argstr[1] == "-")))) pinfo.usage();

    opt.cmdline = ad.cmdline;

    bool is_min_qscore_set(false);
    bool is_min_sascore_set(false);
    bool is_min_pascore_set(false);
    bool is_contigs_set(false);
    bool is_contig_reads_set(false);

    bool is_bsnp_ssd_no_mismatch(false);
    bool is_bsnp_ssd_one_mismatch(false);

    bool is_min_can_indel_set(false);
    bool is_min_can_indel_frac_set(false);
    bool is_min_small_can_indel_frac_set(false);
    bool is_max_can_indel_density_set(false);

    bool is_max_basecall_filter_fraction(false);

    bool is_max_vexp_iterations(false);

    bool is_max_indel_size(false);

    bool is_cif(false);

    bool is_inmp(false);

    bool is_uacmq_set(false);

    const unsigned as(ad.size());
    for (unsigned i(0); i<as; ++i) {
        if (ad.argmark[i]) continue;

        if (ad.argstr[i]=="-lsnp-alpha") {
            set_xrange_arg(i,ad,opt.is_lsnp,opt.lsnp_alpha);
        } else if (ad.argstr[i]=="-bsnp-diploid-file") {
            set_filename_arg(i,ad,opt.is_bsnp_diploid_file,opt.bsnp_diploid_filename);
        } else if (ad.argstr[i]=="-bsnp-diploid-allele-file") {
            set_filename_arg(i,ad,opt.is_bsnp_diploid_allele_file,opt.bsnp_diploid_allele_filename);
        } else if (ad.argstr[i]=="-bsnp-monoploid") {
            set_xrange_arg(i,ad,opt.is_bsnp_monoploid,opt.bsnp_monoploid_theta);
        } else if (ad.argstr[i]=="-bsnp-ssd-no-mismatch") {
            set_xrange_arg(i,ad,is_bsnp_ssd_no_mismatch,opt.bsnp_ssd_no_mismatch,true);
        } else if (ad.argstr[i]=="-bsnp-ssd-one-mismatch") {
            set_xrange_arg(i,ad,is_bsnp_ssd_one_mismatch,opt.bsnp_ssd_one_mismatch,true);
        } else if (ad.argstr[i]=="-bsnp-diploid-het-bias") {
            set_xrange_arg(i,ad,opt.is_bsnp_diploid_het_bias,opt.bsnp_diploid_het_bias,true);
        } else if (ad.argstr[i]=="-bsnp-nploid") {
            set_nploid_arg(i,ad,opt.is_bsnp_nploid,opt.bsnp_nploid_ploidy,opt.bsnp_nploid_snp_prob);
        } else if (ad.argstr[i]=="-bindel-diploid-het-bias") {
            set_xrange_arg(i,ad,opt.is_bindel_diploid_het_bias,opt.bindel_diploid_het_bias,true);
        } else if (ad.argstr[i]=="-bindel-diploid-file") {
            set_filename_arg(i,ad,opt.is_bindel_diploid_file,opt.bindel_diploid_filename);
        } else if (ad.argstr[i]=="-anom-distro-table-alpha") {
            set_xrange_arg(i,ad,opt.is_adis_table,opt.adis_table_alpha);
        } else if (ad.argstr[i]=="-anom-distro-lrt-alpha") {
            set_xrange_arg(i,ad,opt.is_adis_lrt,opt.adis_lrt_alpha);
        } else if (ad.argstr[i]=="-anom-cov-alpha") {
            set_xrange_arg(i,ad,opt.is_acov,opt.acov_alpha);
        } else if (ad.argstr[i]=="-filter-anom-calls") {
            opt.is_filter_anom_calls=true;
        } else if (ad.argstr[i]=="-print-evidence") {
            opt.is_print_evidence=true;
        } else if (ad.argstr[i]=="-print-all-site-evidence") {
            opt.is_print_all_site_evidence=true;
        } else if (ad.argstr[i]=="-min-qscore") {
            set_arg(i,ad,is_min_qscore_set,opt.min_qscore);
        } else if (ad.argstr[i]=="-max-window-mismatch") {
            int max_win_mismatch_tmp;
            set_win_arg(i,ad,opt.is_max_win_mismatch,max_win_mismatch_tmp,opt.max_win_mismatch_flank_size);
            if (max_win_mismatch_tmp<0) {
                pinfo.usage("first argument following -max-window-mismatch must be a non-negative integer\n");
            }
            opt.max_win_mismatch=max_win_mismatch_tmp;
        } else if (ad.argstr[i]=="-min-single-align-score") {
            set_arg(i,ad,is_min_sascore_set,opt.min_single_align_score);
        } else if (ad.argstr[i]=="-min-paired-align-score") {
            set_arg(i,ad,is_min_pascore_set,opt.min_paired_align_score);
        } else if (ad.argstr[i]=="-single-align-score-exclude-mode") {
            opt.single_align_score_exclude_mode=true;
        } else if (ad.argstr[i]=="-single-align-score-rescue-mode") {
            opt.single_align_score_rescue_mode=true;
        } else if (ad.argstr[i]=="-counts") {
            set_filename_arg(i,ad,opt.is_counts,opt.counts_filename);
        } else if (ad.argstr[i]=="-clobber") {
            opt.is_clobber = true;
        } else if (ad.argstr[i]=="-bam-file") {
            bool is_bam_filename(! opt.bam_filename.empty());
            set_filename_arg(i,ad,is_bam_filename,opt.bam_filename);
        } else if (ad.argstr[i]=="-bam-seq-name") {
            bool is_bam_seq_name(! opt.bam_seq_name.empty());
            set_filename_arg(i,ad,is_bam_seq_name,opt.bam_seq_name);
        } else if (ad.argstr[i]=="-samtools-reference") {
            set_filename_arg(i,ad,opt.is_samtools_ref_set,opt.samtools_ref_seq_file);
        } else if (ad.argstr[i]=="-indel-contigs") {
            set_filename_arg(i,ad,is_contigs_set,opt.indel_contig_filename);
        } else if (ad.argstr[i]=="-indel-contig-reads") {
            set_filename_arg(i,ad,is_contig_reads_set,opt.indel_contig_read_filename);
        } else if (ad.argstr[i]=="-indel-error-rate") {
            set_xrange_arg(i,ad,opt.is_simple_indel_error,opt.simple_indel_error,true);
        } else if (ad.argstr[i]=="-indel-nonsite-match-prob") {
            set_xrange_arg(i,ad,is_inmp,opt.indel_nonsite_match_prob,true);
        } else if (ad.argstr[i]=="-report-range-begin") {
            set_arg(i,ad,opt.user_report_range.is_begin_pos,opt.user_report_range.begin_pos);
        } else if (ad.argstr[i]=="-report-range-end") {
            set_arg(i,ad,opt.user_report_range.is_end_pos,opt.user_report_range.end_pos);
        } else if (ad.argstr[i]=="-report-range-reference") {
            opt.is_report_range_ref=true;
        } else if (ad.argstr[i]=="-genome-size") {
            set_arg(i,ad,opt.is_user_genome_size,opt.user_genome_size);
        } else if (ad.argstr[i]=="-print-all-poly-gt") {
            opt.is_print_all_poly_gt=true;
        } else if (ad.argstr[i]=="-print-used-allele-counts") {
            opt.is_print_used_allele_counts=true;
        } else if (ad.argstr[i]=="-used-allele-count-min-qscore") {
            set_arg(i,ad,is_uacmq_set,opt.used_allele_count_min_qscore);
        } else if (ad.argstr[i]=="-min-candidate-indel-reads") {
            set_arg(i,ad,is_min_can_indel_set,opt.default_min_candidate_indel_reads);
        } else if (ad.argstr[i]=="-min-candidate-indel-read-frac") {
            set_xrange_arg(i,ad,is_min_can_indel_frac_set,opt.min_candidate_indel_read_frac,true);
        } else if (ad.argstr[i]=="-min-small-candidate-indel-read-frac") {
            set_xrange_arg(i,ad,is_min_small_can_indel_frac_set,opt.default_min_small_candidate_indel_read_frac,true);
        } else if (ad.argstr[i]=="-max-candidate-indel-density") {
            set_xrange_arg(i,ad,is_max_can_indel_density_set,opt.max_candidate_indel_density,false,true);
        } else if (ad.argstr[i]=="-realigned-read-file") {
            set_filename_arg(i,ad,opt.is_realigned_read_file,opt.realigned_read_filename);
        } else if (ad.argstr[i]=="-realign-submapped-reads") {
            opt.is_realign_submapped_reads=true;
        } else if (ad.argstr[i]=="-snp-max-basecall-filter-fraction") {
            set_xrange_arg(i,ad,is_max_basecall_filter_fraction,opt.max_basecall_filter_fraction,true);
        } else if (ad.argstr[i]=="-max-vexp-iterations") {
            set_arg(i,ad,is_max_vexp_iterations,opt.max_vexp_iterations);
        } else if (ad.argstr[i]=="-min-vexp") {
            set_xrange_arg(i,ad,opt.is_min_vexp,opt.min_vexp,true);
        } else if (ad.argstr[i]=="-no-ambiguous-path-clip") {
            opt.is_clip_ambiguous_path=false;
        } else if (ad.argstr[i]=="-max-indel-size") {
            set_arg(i,ad,is_max_indel_size,opt.max_indel_size);
        } else if (ad.argstr[i]=="-all-warnings") {
            opt.verbosity=LOG_LEVEL::ALLWARN;
        } else if (ad.argstr[i]=="-filter-unanchored") {
            opt.is_filter_unanchored=true;
        } else if (ad.argstr[i]=="-include-singleton") {
            opt.is_include_singleton = true;
        } else if (ad.argstr[i]=="-include-anomalous") {
            opt.is_include_anomalous = true;
        } else if (ad.argstr[i]=="-candidate-indel-file") {
            set_filename_arg(i,ad,is_cif,opt.candidate_indel_filename);
        } else if (ad.argstr[i]=="-write-candidate-indels-only") {
            opt.is_write_candidate_indels_only=true;
        } else if (ad.argstr[i]=="-skip-variable-metadata") {
            opt.is_write_variable_metadata=false;
        } else if (ad.argstr[i]=="-baby-elephant") {
            opt.is_baby_elephant=true;
        } else if (ad.argstr[i]=="-h") {
            pinfo.usage();
        } else { continue; }

        ad.argmark[i] = true;
    }

    ad.finalize_args();


    // although we've tried to separate parameter parsing and
    // validation, this validation code is stuck here for now:
    //

    // sanity check argument settings:
    //
    {
        // sanity check input read info:

        if (opt.bam_filename.empty()) {
            pinfo.usage("Must specify a sorted BAM file containing aligned sample reads");
        }

        if ((! opt.bam_filename.empty()) && opt.bam_seq_name.empty()) {
            pinfo.usage("must specify -bam-seq-name when a bam file is provided");
        }
    }

    if (! opt.is_ref_set()) {
        pinfo.usage("must specify either single-seq-reference or samtools-reference");
    }

    if (opt.is_samtools_ref_set && opt.bam_seq_name.empty()) {
        pinfo.usage("must specify bam-seq-name when using samtools-reference");
    }

    validate_blt_opt(pinfo,opt);
}
