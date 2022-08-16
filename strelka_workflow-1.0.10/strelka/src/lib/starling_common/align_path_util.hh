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

#ifndef __ALIGN_PATH_UTIL_HH
#define __ALIGN_PATH_UTIL_HH

#include "starling_common/align_path.hh"
#include "blt_util/blt_types.hh"
#include "starling_common/starling_types.hh"


#include <cassert>


namespace ALIGNPATH {

inline
void
increment_path(const path_t& path,
               unsigned& path_index,
               unsigned& read_offset,
               pos_t& ref_offset)
{
    const path_segment& ps(path[path_index]);

    if       (ps.type == MATCH) {
        read_offset += ps.length;
        ref_offset += ps.length;
    } else if (ps.type == DELETE || ps.type == SKIP) {
        ref_offset += ps.length;
    } else if (ps.type == INSERT || ps.type == SOFT_CLIP) {
        read_offset += ps.length;
    } else if (ps.type == HARD_CLIP || ps.type == PAD) {
        // do nothing
    } else {
        assert(0); // can't handle other CIGAR types yet
    }

    path_index++;
}


// Initialize to the segment count, insert and delete size of a
// swap in the path. assumes path_index points to the begining of
// a swap:
//
struct swap_info {
    swap_info(const path_t& path,
              const unsigned path_index)
        : n_seg(path_index)
        , insert_length(0)
        , delete_length(0)
    {
        const unsigned aps(path.size());
        for (; (n_seg<aps) && is_segment_type_indel(path[n_seg].type); ++n_seg) {
            const path_segment& ps(path[n_seg]);
            if     (ps.type==INSERT) { insert_length += ps.length; }
            else if (ps.type==DELETE) { delete_length += ps.length; }
            else                     { assert(0); }
        }
        n_seg -= path_index;
    }

    unsigned n_seg;
    unsigned insert_length;
    unsigned delete_length;
};
}


#endif
