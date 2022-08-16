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



#include "starling_pos_processor_indel_util.hh"
#include "starling_read_util.hh"

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "starling_common/align_path_util.hh"

#include "boost/foreach.hpp"

#include <cassert>

#include <iostream>
#include <sstream>



static
void
finish_indel_sppr(indel_observation& obs,
                  starling_pos_processor_base& sppr,
                  const unsigned sample_no) {

    const bool is_novel_indel(sppr.insert_indel(obs,sample_no));

    // contig reads are supposed to be associated with indels from their contig only:
    //
    if ((INDEL_ALIGN_TYPE::CONTIG_READ == obs.data.iat) && is_novel_indel) {
        std::ostringstream oss;
        oss << "ERROR: contig read contains novel indel: " << obs.key << "\n";
        throw blt_exception(oss.str().c_str());
    }
}



static
void
bam_seq_to_str(const bam_seq_base& bs,
               const unsigned start,
               const unsigned end,
               std::string& s) {
    s.clear();
    for (unsigned i(start); i<end; ++i) s.push_back(bs.get_char(i));
}



const unsigned max_cand_filter_insert_size(10);



// Handle edge inserts. For contigs we take the edge insert verbatim. For contig
// reads we need to find out what (usually larger) indel this edge corresponds to.
// for regular reads we record this as a private indel (ie. it does not count towards
// candidacy.
//
static
void
process_edge_insert(const unsigned max_indel_size,
                    const ALIGNPATH::path_t& path,
                    const bam_seq_base& bseq,
                    starling_pos_processor_base& sppr,
                    indel_observation& obs,
                    const unsigned sample_no,
                    const indel_set_t* edge_indel_ptr,
                    const unsigned seq_len,
                    const std::pair<unsigned,unsigned>& ends,
                    const unsigned path_index,
                    const unsigned read_offset,
                    const pos_t ref_head_pos,
                    const bool is_pinned_indel)
{
    using namespace ALIGNPATH;

    const path_segment& ps(path[path_index]);

    if       (obs.data.iat == INDEL_ALIGN_TYPE::CONTIG) {
        obs.key.pos=ref_head_pos;
        obs.key.length=ps.length;
        if (path_index==ends.first) { // right side BP:
            obs.key.type=INDEL::BP_RIGHT;
            const unsigned next_read_offset(read_offset+ps.length);
            const unsigned start_offset(next_read_offset-std::min(next_read_offset,max_indel_size));
            bam_seq_to_str(bseq,start_offset,next_read_offset,obs.data.insert_seq);
        } else { // left side BP:
            obs.key.type=INDEL::BP_LEFT;
            const unsigned seq_size(std::min(seq_len-read_offset,max_indel_size));
            bam_seq_to_str(bseq,read_offset,read_offset+seq_size,obs.data.insert_seq);
        }
        finish_indel_sppr(obs,sppr,sample_no);

    } else if (obs.data.iat == INDEL_ALIGN_TYPE::CONTIG_READ) {

        // add edge indels for contig reads:
        if (NULL != edge_indel_ptr) {
            const pos_t current_pos(ref_head_pos);
            BOOST_FOREACH(const indel_key& ik, *edge_indel_ptr) {
                // This checks that we've identified the edge indel
                // for this read out of the full set of indels in the
                // contig. It is not fantastically robust. In summary:
                //
                // Contig indel 'j' is identified as our edge indel if
                // the read insertion starts at the same place as the
                // contig indel unless this is the first segment in
                // the read, in which case we check that the read
                // insertion ends at the same place as the contig
                // indel.
                //
                if (path_index!=ends.first) {
                    if (current_pos!=ik.pos) continue;
                } else {
                    if (current_pos!=ik.right_pos()) continue;
                }

                obs.key = ik;

                // large insertion breakpoints are not filtered as noise:
                if (obs.data.is_noise) {
                    if (obs.key.is_breakpoint() ||
                        ((obs.key.type == INDEL::INSERT) &&
                         (obs.key.length > max_cand_filter_insert_size))) {
                        obs.data.is_noise=false;
                    }
                }

                finish_indel_sppr(obs,sppr,sample_no);
                break;
            }
        }
    } else {
        // do not allow edge-indels on genomic reads to generate or support candidate
        // indels, except for pinned cases:
        if (! is_pinned_indel) return;
        if (ps.length > max_indel_size) return;

        obs.key.pos=ref_head_pos;
        obs.key.length = ps.length;
        assert(ps.type == INSERT);
        obs.key.type=INDEL::INSERT;

#ifdef SPI_DEBUG
        log_os << "FOOBAR: adding pinned edge indel: " << obs.key << "\n";
#endif

        bam_seq_to_str(bseq,read_offset,read_offset+ps.length,obs.data.insert_seq);
        finish_indel_sppr(obs,sppr,sample_no);
    }
}



/// handle edge deletions
///
/// currently these are only alllowed for indels adjacent to an intron
///
static
void
process_edge_delete(const unsigned max_indel_size,
                    const ALIGNPATH::path_t& path,
                    const bam_seq_base& bseq,
                    starling_pos_processor_base& sppr,
                    indel_observation& obs,
                    const unsigned sample_no,
                    const unsigned path_index,
                    const unsigned read_offset,
                    const pos_t ref_head_pos,
                    const bool is_pinned_indel)
{
    // we should never see grouper reads here
    if ((obs.data.iat == INDEL_ALIGN_TYPE::CONTIG) ||
        (obs.data.iat == INDEL_ALIGN_TYPE::CONTIG_READ)) {
        return;
    }


    using namespace ALIGNPATH;

    const path_segment& ps(path[path_index]);

    // do not allow edge-indels on genomic reads to generate or support candidate
    // indels, except for pinned cases:
    if (! is_pinned_indel) return;
    if (ps.length > max_indel_size) return;

    obs.key.pos=ref_head_pos;
    obs.key.length = ps.length;
    assert(ps.type == DELETE);
    obs.key.type=INDEL::DELETE;

#ifdef SPI_DEBUG
    log_os << "FOOBAR: adding pinned edge indel: " << obs.key << "\n";
#endif

    bam_seq_to_str(bseq,read_offset,read_offset+ps.length,obs.data.insert_seq);
    finish_indel_sppr(obs,sppr,sample_no);
}



// Note that unlike other indel processors, the swap processor returns
// the number of segments consumed.
//
// Like regular indels, a swap will be reported as a breakpoint if
// it's too long (meaning longeset of insert,delete size).
//
static
unsigned
process_swap(const unsigned max_indel_size,
             const ALIGNPATH::path_t& path,
             const bam_seq_base& bseq,
             starling_pos_processor_base& sppr,
             indel_observation& obs,
             const unsigned sample_no,
             const unsigned path_index,
             const unsigned read_offset,
             const pos_t ref_head_pos)
{
    using namespace ALIGNPATH;

    const swap_info sinfo(path,path_index);
    const unsigned swap_size(std::max(sinfo.insert_length,sinfo.delete_length));

    // large insertions are not filtered as noise:
    if (obs.data.is_noise) {
        if (sinfo.insert_length > max_cand_filter_insert_size) {
            obs.data.is_noise=false;
        }
    }

    if (swap_size <= max_indel_size) {
        obs.key.pos=ref_head_pos;
        obs.key.length=sinfo.insert_length;
        obs.key.swap_dlength=sinfo.delete_length;
        obs.key.type = INDEL::SWAP;
        bam_seq_to_str(bseq,read_offset,read_offset+sinfo.insert_length,obs.data.insert_seq);
        finish_indel_sppr(obs,sppr,sample_no);

    } else {

        // left side BP:
        {
            obs.key.pos=ref_head_pos;
            obs.key.length=swap_size;
            obs.key.type=INDEL::BP_LEFT;
            const unsigned start(read_offset);
            const unsigned size(bseq.size()-read_offset);
            const unsigned end(start+std::min(size,max_indel_size));
            bam_seq_to_str(bseq,start,end,obs.data.insert_seq);
            finish_indel_sppr(obs,sppr,sample_no);
        }

        // right side BP:
        {
            obs.key.pos=ref_head_pos+sinfo.delete_length;
            obs.key.length=swap_size;
            obs.key.type=INDEL::BP_RIGHT;
            const unsigned next_read_offset(read_offset+sinfo.insert_length);
            const unsigned start_offset(next_read_offset-std::min(next_read_offset,max_indel_size));
            bam_seq_to_str(bseq,start_offset,next_read_offset,obs.data.insert_seq);
            finish_indel_sppr(obs,sppr,sample_no);
        }
    }

    return sinfo.n_seg;
}



// Handle the regular ol' insertions and deletions. Reports these
// types as breakpoints when they're too long:
//
static
void
process_simple_indel(const unsigned max_indel_size,
                     const ALIGNPATH::path_t& path,
                     const bam_seq_base& bseq,
                     starling_pos_processor_base& sppr,
                     indel_observation& obs,
                     const unsigned sample_no,
                     const unsigned path_index,
                     const unsigned read_offset,
                     const pos_t ref_head_pos)
{
    using namespace ALIGNPATH;

    const path_segment& ps(path[path_index]);

    // large insertion breakpoints are not filtered as noise:
    if (obs.data.is_noise) {
        if ((ps.type == INSERT) &&
            (ps.length > max_cand_filter_insert_size)) {
            obs.data.is_noise=false;
        }
    }

    if (ps.length <= max_indel_size) {
        obs.key.pos=ref_head_pos;
        obs.key.length = ps.length;
        if (ps.type == INSERT) {
            obs.key.type=INDEL::INSERT;
            bam_seq_to_str(bseq,read_offset,read_offset+ps.length,obs.data.insert_seq);
        } else {
            obs.key.type=INDEL::DELETE;
        }
        finish_indel_sppr(obs,sppr,sample_no);
    } else {
        // left side BP:
        {
            obs.key.pos=ref_head_pos;
            obs.key.length=ps.length;
            obs.key.type=INDEL::BP_LEFT;
            const unsigned start(read_offset);
            const unsigned size(bseq.size()-read_offset);
            const unsigned end(start+std::min(size,max_indel_size));
            bam_seq_to_str(bseq,start,end,obs.data.insert_seq);
            finish_indel_sppr(obs,sppr,sample_no);
        }
        // right side BP:
        {
            obs.key.pos=ref_head_pos;
            if (ps.type == DELETE) obs.key.pos+=ps.length;
            obs.key.length=ps.length;
            obs.key.type=INDEL::BP_RIGHT;

            const unsigned next_read_offset(read_offset+((ps.type==INSERT) ? ps.length : 0));
            const unsigned start_offset(next_read_offset-std::min(next_read_offset,max_indel_size));
            bam_seq_to_str(bseq,start_offset,next_read_offset,obs.data.insert_seq);
            finish_indel_sppr(obs,sppr,sample_no);
        }
    }
}



// Extract indel information from an alignment and store this
// in the starling_pos_processor indel buffer.
//
// assumes that path is already validated for seq!!!
//
unsigned
add_alignment_indels_to_sppr(const unsigned max_indel_size,
                             const reference_contig_segment& ref,
                             const alignment& al,
                             const bam_seq_base& read_seq,
                             starling_pos_processor_base& sppr,
                             const INDEL_ALIGN_TYPE::index_t iat,
                             const align_id_t id,
                             const unsigned sample_no,
                             const std::pair<bool,bool>& edge_pin,
                             const indel_set_t* edge_indel_ptr) {

    using namespace ALIGNPATH;

    const unsigned seq_len(read_seq.size());

    if (is_apath_invalid(al.path,seq_len)) {
        std::ostringstream oss;
        oss << "ERROR: Can't handle alignment path '" << apath_to_cigar(al.path) << "' -- " << get_apath_invalid_reason(al.path,seq_len) << "\n";
        throw blt_exception(oss.str().c_str());
    }

    if (is_apath_starling_invalid(al.path)) {
        std::ostringstream oss;
        oss << "ERROR: can't handle alignment path '" << apath_to_cigar(al.path) << "'\n";
        throw blt_exception(oss.str().c_str());
    }

    const rc_segment_bam_seq ref_bseq(ref);

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(al.path));

    pos_range valid_pr;
    get_valid_alignment_range(al,ref_bseq,read_seq,valid_pr);

    unsigned path_index(0);
    unsigned read_offset(0);
    pos_t ref_head_pos(al.pos);

    unsigned total_indel_ref_span_per_read(0);

    const unsigned aps(al.path.size());
    while (path_index<aps) {
        const path_segment& ps(al.path[path_index]);
        const bool is_begin_edge(path_index<ends.first);
        const bool is_end_edge(path_index>ends.second);
        const bool is_edge_segment(is_begin_edge || is_end_edge);

        const bool is_swap_start(is_segment_swap_start(al.path,path_index));

        assert(ps.type != SKIP);
        assert(! (is_edge_segment && is_swap_start));

        indel_observation obs;
        obs.data.iat = iat;
        obs.data.id = id;

        if (MATCH != ps.type) {
            pos_range indel_read_pr;
            indel_read_pr.set_begin_pos((read_offset==0) ? 0 : (read_offset-1));

            unsigned rlen(0);
            if       (is_swap_start) {
                const swap_info sinfo(al.path,path_index);
                rlen=sinfo.insert_length;

                if (sinfo.delete_length<=max_indel_size) {
                    total_indel_ref_span_per_read += sinfo.delete_length;
                }
            } else if (is_segment_type_read_length(ps.type)) {
                rlen=ps.length;
            } else {
                if (ps.type == DELETE) {
                    if (ps.length <= max_indel_size) {
                        total_indel_ref_span_per_read += ps.length;
                    }
                }
            }
            indel_read_pr.set_end_pos(std::min(seq_len,read_offset+1+rlen));
            if (! valid_pr.is_superset_of(indel_read_pr)) obs.data.is_noise=true;
        }

        unsigned n_seg(1); // number of path segments consumed
        if (is_edge_segment) {
            // is this indel occurring on a pinned edge (ie against an exon?)
            const bool is_pinned_indel((is_begin_edge && edge_pin.first) ||
                                       (is_end_edge && edge_pin.second));

            // edge inserts are allowed for intron adjacent and grouper reads, edge deletions for intron adjacent only:
            if (ps.type == INSERT) {
                process_edge_insert(max_indel_size,al.path,read_seq,
                                    sppr,obs,sample_no,edge_indel_ptr,
                                    seq_len,ends,path_index,read_offset,ref_head_pos,
                                    is_pinned_indel);
            } else if (ps.type == DELETE) {
                if (is_pinned_indel) {
                    process_edge_delete(max_indel_size,al.path,read_seq,
                                        sppr,obs,sample_no,
                                        path_index,read_offset,ref_head_pos,
                                        is_pinned_indel);
                }
            }
        } else if (is_swap_start) {
            n_seg = process_swap(max_indel_size,al.path,read_seq,
                                 sppr,obs,sample_no,
                                 path_index,read_offset,ref_head_pos);

        } else if (is_segment_type_indel(al.path[path_index].type)) {
//            log_os << "offset " << read_offset <<  endl;
//            log_os << " " << ref_head_pos <<  endl;
            process_simple_indel(max_indel_size,al.path,read_seq,
                                 sppr,obs,sample_no,
                                 path_index,read_offset,ref_head_pos);

        }

        for (unsigned i(0); i<n_seg; ++i) { increment_path(al.path,path_index,read_offset,ref_head_pos); }
    }

    return total_indel_ref_span_per_read;
}


