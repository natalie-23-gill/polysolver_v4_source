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

#pragma once

#include "blt_util/pos_range.hh"

#include <ciso646>
#include <iosfwd>
#include <stdint.h>
#include <string>
#include <vector>

/// process export MD and output the alignment information in a format
/// which will not need to change significantly for SAM/BAM (by
/// providing essentially a parsed forward strand CIGAR string)
///
/// idea is for this code to migrate to some kind of joint export/sam api code
///
namespace ALIGNPATH {

enum align_t {
    NONE,
    MATCH,
    INSERT,
    DELETE,
    SKIP,
    SOFT_CLIP,
    HARD_CLIP,
    PAD
};

inline
char
segment_type_to_cigar_code(const align_t id) {
    switch (id) {
    case MATCH     : return 'M';
    case INSERT    : return 'I';
    case DELETE    : return 'D';
    case SKIP      : return 'N';
    case SOFT_CLIP : return 'S';
    case HARD_CLIP : return 'H';
    case PAD       : return 'P';
    default :        return 'X';
    }
}

inline
align_t
cigar_code_to_segment_type(const char c) {
    switch (c) {
    case 'M' : return MATCH;
    case 'I' : return INSERT;
    case 'D' : return DELETE;
    case 'N' : return SKIP;
    case 'S' : return SOFT_CLIP;
    case 'H' : return HARD_CLIP;
    case 'P' : return PAD;
    default  : return NONE;
    }
}

inline
bool
is_segment_type_read_length(const align_t id) {
    switch (id) {
    case MATCH     :
    case INSERT    :
    case SOFT_CLIP : return true;
    default        : return false;
    }
}

inline
bool
is_segment_type_ref_length(const align_t id) {
    switch (id) {
    case MATCH  :
    case DELETE :
    case SKIP   : return true;
    default     : return false;
    }
}

inline
bool
is_segment_type_indel(const align_t id) {
    switch (id) {
    case INSERT :
    case DELETE : return true;
    default     : return false;
    }
}

struct path_segment {
    path_segment(const align_t t = NONE,
                 const unsigned l = 0) : type(t), length(l) {}

    void
    clear() {
        type=NONE;
        length=0;
    }

    bool
    operator==(const path_segment& rhs) const {
        return ((type==rhs.type) and (length==rhs.length));
    }

    // arbitrary ordering which lets us look up from a set of alignments:
    bool
    operator<(const path_segment& rhs) const {
        if (type<rhs.type) return true;
        if (type==rhs.type) {
            return length<rhs.length;
        }
        return false;
    }

    align_t type;
    unsigned length;
};

typedef std::vector<path_segment> path_t;

std::ostream& operator<<(std::ostream& os, const path_t& apath);

void
export_md_to_apath(const char* md,
                   const bool is_fwd_strand,
                   path_t& apath,
                   const bool is_edge_deletion_error=true);

void
apath_to_cigar(const path_t& apath,
               std::string& cigar);

inline
std::string
apath_to_cigar(const path_t& apath) {
    std::string cigar;
    apath_to_cigar(apath,cigar);
    return cigar;
}

void
apath_to_export_md(path_t& apath,
                   const char* ref_seq,
                   const char* ref_end,
                   const int32_t ref_pos,
                   const std::string& read_bases,
                   const bool is_fwd_strand,
                   std::string& md);

// Convert cigar string into apath format
//
// any padding in the CIGAR string is removed
void
cigar_to_apath(const char* cigar,
               path_t& apath);

unsigned
apath_read_length(const path_t& apath);

unsigned
apath_ref_length(const path_t& apath);

// how much unaligned sequence (soft_clip or insert) occurs before the first aligned base?
unsigned
apath_read_lead_size(const path_t& apath);

// how much unaligned sequence (soft_clip or insert) occurs after the last aligned base?
unsigned
apath_read_trail_size(const path_t& apath);

// how much soft_clip occurs before the first aligned base?
unsigned
apath_soft_clip_lead_size(const path_t& apath);

// how much soft_clip occurs after the last aligned base?
unsigned
apath_soft_clip_trail_size(const path_t& apath);

// how much insert occurs before the first aligned base?
unsigned
apath_insert_lead_size(const path_t& apath);

// how much insert occurs after the last aligned base?
unsigned
apath_insert_trail_size(const path_t& apath);

// remove any edge clip from apath and return the amount
// removed from each side. if ambiguous, lead is favored over trail
void
apath_clip_clipper(path_t& apath,
                   unsigned& hc_lead,
                   unsigned& hc_trail,
                   unsigned& sc_lead,
                   unsigned& sc_trail);

// adds lead clip to front of alignment and trail clip
// to back -- assumes no clipping exists on the path already.
//
void
apath_clip_adder(path_t& apath,
                 const unsigned hc_lead,
                 const unsigned hc_trail,
                 const unsigned sc_lead,
                 const unsigned sc_trail);

// 'cleans' the path so that it can be used, or used more consistently, in starka
//
// note this does not try to correct or work around anything
// which can't be unambiguously reinterpreted to a simpler form
//
// 1. remove zero length alignment segments
// 2. remove pad segments
// 3. remove repeated segments
// 4. for any combined insertion/deletion pair, reduce this to
//    a single segment pair (in either order)
//
//  \return true if path has been altered
//
bool
apath_cleaner(path_t& apath);

#if 0
// Get the match descriptor segment numbers for the first and last
// non-soft/hard clipped segments. Return total segment size on
// error.
std::pair<unsigned,unsigned>
get_nonclip_end_segments(const path_t& apath);
#endif

// return the read coordinate range after clipping:
pos_range
get_nonclip_range(const path_t& apath);

// Get the match descriptor segment numbers for the first and last
// match segments. Return total segment size on error.
std::pair<unsigned,unsigned>
get_match_edge_segments(const path_t& apath);

unsigned
apath_exon_count(const path_t& apath);

// provide reference offsets for the begining of each exon:
//
struct exon_offsets {
    exon_offsets(const path_t& apath)
        : _apath(apath)
        , _asize(apath.size())
        , _offset(0)
        , _segment(0)
    {}

    bool
    next() {
        bool is_break_next(false);
        for (; _segment<_asize; ++_segment) {
            if (is_break_next) return true;
            const path_segment& ps(_apath[_segment]);
            if (ps.type==SKIP) is_break_next=true;
            if (is_segment_type_ref_length(ps.type)) _offset += ps.length;
        }
        return false;
    }

    unsigned
    offset() const { return _offset; }

private:
    const path_t& _apath;
    const unsigned _asize;
    unsigned _offset;
    unsigned _segment;
};

/// does the alignment contain any soft-clipped segments?
bool
is_soft_clipped(const path_t& apath);

/// is either edge of the alignment soft-clipped or hard-clipped?
bool
is_clipped(const path_t& apath);

/// does either edge of the alignment
/// contain a segment which impacts read length or reference positions?
/// (INSERT,DELETE,SKIP,SOFT_CLIP)
///
/// Note: "edge" is defined as any segment with match segments to only one side
/// Note: edge HARD_CLIP, PAD, etc.. are ignored
///
bool
is_edge_readref_len_segment(const path_t& apath);

// does alignment contain an adjacent insertion/deletion event?
//
bool
is_seq_swap(const path_t& apath);

// is the given segment the begining of a seq swap?
bool
is_segment_swap_start(const path_t& apath,
                      const unsigned i);

// test if alignment has no match:
bool
is_apath_floating(const path_t& apath);


namespace ALIGN_ISSUE {
enum issue_t {
    NONE,
    CLIPPING,
    EDGE_DELETE,
    EDGE_SKIP,
    UNKNOWN_SEGMENT,
    REPEATED_SEGMENT,
    FLOATING,
    LENGTH
};

inline
const char*
description(const issue_t i) {
    switch (i) {
    case CLIPPING: return "alignment contains invalid clipping";
    case EDGE_DELETE: return "deletion on alignment edge";
    case EDGE_SKIP: return "skip on alignment edge";
    case UNKNOWN_SEGMENT: return "unknown segment in alignment";
    case REPEATED_SEGMENT: return "alignment contains repeated segment";
    case FLOATING: return "alignment contains no match segments";
    case LENGTH: return "alignment length does not match read length";
    default: return "no error";
    }
}
}


/// Take a shot at the relatively simple stuff:
///
/// 1) clipping only occurs on the edge and hardclip must occur outside of soft-clip
/// 2) delete and skip cannot occur on edge
///   2a) delete and skip cannot occur with only insert and clip connecting them to edge
/// 3) no unknown segments
/// 4) no repeated segments
///      Note this might semi-legitimately occur where padding is stripped out of an alignment.
/// 5) must contain at least one match segment
///
ALIGN_ISSUE::issue_t
get_apath_invalid_type(const path_t& path,
                       const unsigned seq_length);

/// if is_apath_invalid fails, this supplies an error string
std::string
get_apath_invalid_reason(const path_t& apath,
                         const unsigned seq_length);

/// simple boolean call to the invalid alignment typer.
inline
bool
is_apath_invalid(const path_t& apath,
                 const unsigned seq_length) {

    return (ALIGN_ISSUE::NONE != get_apath_invalid_type(apath,seq_length));
}

// check for conditions on an otherwise valid path which starling
// does not handle:
bool
is_apath_starling_invalid(const path_t& apath);

#if 0
normalize_path();
#endif
}

