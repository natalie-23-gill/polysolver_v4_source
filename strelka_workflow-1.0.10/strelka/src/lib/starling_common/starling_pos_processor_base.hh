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

#pragma once


#include "blt_common/adjust_joint_eprob.hh"
#include "blt_common/map_level.hh"
#include "blt_util/depth_stream_stat_range.hh"
#include "blt_util/pos_processor_base.hh"
#include "blt_util/stage_manager.hh"
#include "blt_util/window_util.hh"
#include "starling_common/depth_buffer.hh"
#include "starling_common/gvcf_aggregator.hh"
#include "starling_common/indel_buffer.hh"
#include "starling_common/indel_set.hh"
#include "starling_common/indel_synchronizer.hh"
#include "starling_common/pos_basecall_buffer.hh"
#include "starling_common/read_mismatch_info.hh"
#include "starling_common/starling_read_buffer.hh"
#include "starling_common/starling_shared.hh"
#include "starling_common/starling_streams_base.hh"


#include "boost/utility.hpp"

#include <memory>
#include <string>

struct diploid_genotype;
struct nploid_info;


//int
//get_influence_zone_size(const unsigned max_indel_size);


/// keep a single copy of this struct to reuse for every site to lower alloc costs:
struct extra_position_data {
    /// stores the column of basecalls actually used for snp-calling after the
    /// mismatch density filter and other quality filters have been applied:
    snp_pos_info good_pi;

    /// stores information on approximate qscore reductions implemented to represent
    /// site-specific basecalling dependency, note this is not applied in somatic
    /// calling:
    std::vector<float> dependent_eprob;
};


/// \brief accumulate sequential position specific information and
/// send to a snp-calling/indel-calling routine after all position
/// information is found
///
/// pos_processor assumes that information related to each position
/// will be available in an approximately sequential fashion, where all
/// position values submitted after position X will be greater than
/// X-POS_BUFFER_SIZE+1. A violation of this assumption will trigger a
/// runtime error.
///
/// The implementation should be split into a generic semi-sequential
/// position-processor and an object with the application-specific
/// position code (here, snp-calling). Maybe once a second application
/// comes up...
///
///
/// range policy:
///
/// if begin_pos is not specified, then event processiing and
/// reporting start at the first pos >= 0 with snp/indel information
/// submitted, else at begin_pos
///
/// if end_pos is not specified, then event processing and reporting
/// end after last_pos with snp/indel information submitted, else at
/// end_pos. If the reference sequence is set then end_pos must be <
/// ref_size and will be adjusted accordingly.
///
/// Submission of snps < first_pos will be ignored. Submission of indels
/// between first_pos and first_pos-MAX_READ_SIZE will be processed but
/// not reported.
///
/// ...
///
struct starling_pos_processor_base : public pos_processor_base, private boost::noncopyable {

    typedef pos_processor_base base_t;

    starling_pos_processor_base(const starling_options& client_opt,
                                const starling_deriv_options& client_dopt,
                                const reference_contig_segment& ref,
                                const starling_streams_base& client_io,
                                const unsigned n_samples);

    virtual
    ~starling_pos_processor_base();

    // finish position report and reset structure to ground state:
    //
    void reset();

    // note that indel position should be normalized before calling:
    //
    // returns true if this indel is novel to the buffer
    //
    bool
    insert_indel(const indel_observation& obs,
                 const unsigned sample_no);

    unsigned
    get_estimated_depth(const pos_t pos,
                        const unsigned sample_no) const;


    // in range [begin,end), is the estimated depth always below
    // depth?
    bool
    is_estimated_depth_range_ge_than(const pos_t begin,
                                     const pos_t end,
                                     const unsigned depth,
                                     const unsigned sample_no) const;

    // first return value is true if the alignment is accepted into
    // the buffer (alignments can fail a number of quality checks --
    // such as being located too far away from other alignments of the
    // same read or having an indel that is too large
    //
    // second return value is read_id
    //
    std::pair<bool,align_id_t>
    insert_read(const bam_record& br,
                const alignment& al,
                const READ_ALIGN::index_t rat,
                const char* chrom_name,
                const MAPLEVEL::index_t maplev,
                const unsigned sample_no,
                const align_id_t contig_id = 0,
                const indel_set_t* contig_indels_ptr = NULL);

#if 0
    starling_read*
    get_read(const align_id_t read_id,
             const unsigned sample_no) {
        return sample(sample_no).read_buff.get_read(read_id);
    }
#endif

    void
    set_head_pos(const pos_t pos);

    // Is a read potentially containing an indel or an indel itself
    // far enough from the report start and stop positions to be
    // excluded?
    //
    bool
    is_range_outside_report_influence_zone(const pos_range& pr) const {
        return (! _report_influence_range.is_range_intersect(pr));
    }

    // Does a range fall outside of the report start and stop positions?
    //
    bool
    is_range_outside_report_zone(const pos_range& pr) const {
        return (! _client_dopt.report_range_limit.is_range_intersect(pr));
    }

protected:
    std::ostream*
    get_report_osptr() const {
        return _client_io.report_osptr();
    }


    void
    process_pos_snp_single_sample(const pos_t pos,
                                  const unsigned sample_no);

    void
    process_pos_indel_single_sample(const pos_t pos,
                                    const unsigned sample_no);

    /////////////////////////////////
    struct win_avg_set {

        win_avg_set(const unsigned size)
            : ss_used_win(size)
            , ss_filt_win(size)
            , ss_spandel_win(size)
            , ss_submap_win(size)
        {}

        window_average ss_used_win;
        window_average ss_filt_win;
        window_average ss_spandel_win;
        window_average ss_submap_win;
    };


    struct win_avgs {

        win_avgs(const starling_options& opt)
            : _max_winsize(0)
            , _is_last_pos(false)
            , _last_insert_pos(false) {
            const unsigned vs(opt.variant_windows.size());
            _wav.resize(vs);
            for (unsigned i(0); i<vs; ++i) {
                const unsigned winsize(1+opt.variant_windows[i].flank_size*2);
                _wav[i].reset(new win_avg_set(winsize));
                if (winsize>_max_winsize) _max_winsize=winsize;
            }
        }

        void
        insert(const pos_t pos,
               const unsigned n_used,
               const unsigned n_filt,
               const unsigned n_spandel,
               const unsigned n_submap) {
            if (_wav.empty()) return;
            check_skipped_pos(pos);
            insert_impl(n_used,n_filt,n_spandel,n_submap);
        }

        // TODO: does check_skipped_pos need to be called here as well?
        void
        insert_null(const pos_t pos) {
            check_skipped_pos(pos);
            const unsigned ws(_wav.size());
            for (unsigned i(0); i<ws; ++i) {
                _wav[i]->ss_used_win.insert_null();
                _wav[i]->ss_filt_win.insert_null();
                _wav[i]->ss_spandel_win.insert_null();
                _wav[i]->ss_submap_win.insert_null();
            }
        }

        win_avg_set&
        get_win_avg_set(const unsigned i) { return *(_wav[i]); }

    private:

        void
        check_skipped_pos(const pos_t pos) {
            if (_is_last_pos && (pos>(_last_insert_pos+1))) {
                const unsigned rep(std::min(static_cast<pos_t>(_max_winsize),(pos-(_last_insert_pos+1))));
                for (unsigned i(0); i<rep; ++i) insert_impl(0,0,0,0);
            }
            _last_insert_pos=pos;
            _is_last_pos=true;
        }

        void
        insert_impl(const unsigned n_used,
                    const unsigned n_filt,
                    const unsigned n_spandel,
                    const unsigned n_submap) {
            const unsigned ws(_wav.size());
            for (unsigned i(0); i<ws; ++i) {
                _wav[i]->ss_used_win.insert(n_used);
                _wav[i]->ss_filt_win.insert(n_filt);
                _wav[i]->ss_spandel_win.insert(n_spandel);
                _wav[i]->ss_submap_win.insert(n_submap);
            }
        }

        std::vector<boost::shared_ptr<win_avg_set> > _wav;
        unsigned _max_winsize;
        bool _is_last_pos;
        pos_t _last_insert_pos;
    };


    // this structure contains all information which is sample dependent:
    //
    struct sample_info {

        sample_info(const starling_options& opt,
                    const unsigned report_size,
                    const unsigned knownref_report_size,
                    read_id_counter* ricp)
            : indel_buff(opt.max_indel_size)
            , read_buff(ricp)
            , sample_opt(opt)
            , isync_default(indel_buff,estdepth_buff,sample_opt)
            , indel_sync_ptr(&isync_default)
            , ss(report_size)
            , used_ss(report_size)
            , ssn(knownref_report_size)
            , used_ssn(knownref_report_size)
            , wav(opt)
        {}

        indel_synchronizer&
        indel_sync() { return *indel_sync_ptr; }


        indel_buffer indel_buff;
        pos_basecall_buffer bc_buff;
        starling_read_buffer read_buff;
        depth_buffer estdepth_buff; // provide an early estimate of read depth before realignment.

        starling_sample_options sample_opt;

        indel_synchronizer isync_default;
        indel_synchronizer* indel_sync_ptr;

        depth_stream_stat_range ss;
        depth_stream_stat_range used_ss;
        depth_stream_stat_range ssn;
        depth_stream_stat_range used_ssn;

        // regional basecall average windows:
        win_avgs wav;

        // keep a single copy of this struct to reuse for every site to lower alloc costs:
        extra_position_data epd;
    };

    sample_info&
    sample(const unsigned sample_no = 0) { return *_sample[sample_no]; }

    const sample_info&
    sample(const unsigned sample_no = 0) const { return *_sample[sample_no]; }

    bool
    empty() const {
        if (! _is_skip_process_pos) {
            for (unsigned s(0); s<_n_samples; ++s) {
                const sample_info& sif(sample(s));
                if (! sif.read_buff.empty()) return false;
                if (! sif.bc_buff.empty()) return false;
            }
            if (! _variant_print_pos.empty()) return false;
            _is_skip_process_pos=true;
        }
        return true;
    }

    // data for haplotoype regions, shared between all samples:
    //
    struct htype_region_data {
        htype_region_data()
            : is_first_region(true)
            , region_alignment(0) {}

        bool is_first_region;
        pos_t region_alignment;
    };

    ///////////////////////////////
    // static methods:
    //
    static
    void
    report_stream_stat(const depth_stream_stat_range& ss,
                       const char* label,
                       const pos_range& pr,
                       std::ostream& os);

private:
    bool
    is_pos_reportable(const pos_t pos) {
        return _client_dopt.report_range_limit.is_pos_intersect(pos);
    }

    void
    insert_pos_submap_count(const pos_t pos,
                            const unsigned sample_no);

    void
    insert_pos_spandel_count(const pos_t pos,
                             const unsigned sample_no);
    void
    insert_mapq_count(const pos_t pos,
                      const unsigned sample_no,
                      const uint8_t mapq);
    void
    update_ranksum(const pos_t pos,
                   const unsigned sample_no,
                   const base_call& bc,
                   const uint8_t mapq,
                   const int cycle);

    void
    insert_pos_basecall(const pos_t pos,
                        const unsigned sample_no,
                        const bool is_tier1,
                        const base_call& bc);

    void
    insert_hap_cand(const pos_t pos,
                    const unsigned sample_no,
                    const bool is_tier1,
                    const bam_seq_base& read_seq,
                    const uint8_t* qual,
                    const unsigned offset);

    void
    process_pos(const int stage_no,
                const pos_t pos);

    void
    load_read_in_depth_buffer(const read_segment& rseg,
                              const unsigned sample_no);

    void
    init_read_segment(const read_segment& rseg,
                      const unsigned sample_no);

    void
    init_read_segment_pos(const pos_t pos);

    void
    align_pos(const pos_t pos);

    void
    rebuffer_pos_reads(const pos_t pos);

    void
    pileup_pos_reads(const pos_t pos);

    void
    pileup_read_segment(const read_segment& rseg,
                        const unsigned sample_no);

    void
    get_region_haplotypes(const known_pos_range full_pr,
                          const known_pos_range active_pr);

    void
    process_htype_pos(const pos_t begin_pos);

    void
    write_reads(const pos_t pos);

    void
    write_candidate_indels_pos(const pos_t pos);

    void
    process_pos_snp_single_sample_impl(const pos_t pos,
                                       const unsigned sample_no);

    const diploid_genotype&
    get_empty_dgt(const char ref) const;

    //////
    virtual
    void
    process_pos_variants(const pos_t pos) = 0;

    //////
    void
    print_delayed_results(const int stage_no,
                          const pos_t pos);

    virtual
    void
    write_counts(const pos_range& output_report_range) const = 0;

    unsigned
    get_largest_read_size() const {
        return _rmi.size();
    }

    // return false if read is too large
    bool
    update_largest_read_size(const unsigned rs);

    unsigned
    get_largest_total_indel_ref_span_per_read() const {
        return _largest_total_indel_ref_span_per_read;
    }

    void
    update_largest_indel_ref_span(const unsigned is);

    void
    update_largest_total_indel_ref_span_per_read(const unsigned is);

    void
    update_stageman();

protected:

    //////////////////////////////////
    // data:
    //
    const starling_options& _client_opt;
    const starling_deriv_options& _client_dopt;
    const reference_contig_segment& _ref;
    const starling_streams_base& _client_io;

    // read-length data structure used to compute mismatch density filter:
    read_mismatch_info _rmi;

    // largest delete length observed for any one indel (but not greater than max_delete_size)
    unsigned _largest_indel_ref_span;

    // largest
    unsigned _largest_total_indel_ref_span_per_read;

    stage_manager _stageman;

    pos_range _report_influence_range;
    const std::string& _chrom_name;

    // used to keep read id's unique across multiple samples:
    read_id_counter _ric;

    unsigned _n_samples;
    sample_info* _sample[MAX_SAMPLE];

    bool _is_dependent_eprob;
    dependent_prob_cache _dpcache;

    std::auto_ptr<diploid_genotype> _empty_dgt[N_BASE];
    std::auto_ptr<nploid_info> _ninfo;
    double* _ws;

    std::set<pos_t> _variant_print_pos;

    htype_region_data _hregion;

    bool _is_variant_windows;

    std::auto_ptr<gvcf_aggregator> _gvcfer;

    // a caching term used for gvcf:
    site_info _site_info;
};
