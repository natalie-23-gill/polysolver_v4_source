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


#include "starling_read_align_clipper.hh"

#include "blt_util/blt_exception.hh"

#include <cassert>

#include <sstream>


//#define DEBUG_ALIGN_CLIP // extra info on ambiguous path clipping...


#ifdef DEBUG_ALIGN_CLIP
#include <iostream>
#endif


namespace {

struct ref_map_type {
    enum map_t {
        NONE,
        MATCH,
        INSERT,
        SOFT_CLIP,
        CONFLICT
    };


    ref_map_type(const map_t t = NONE,
                 const pos_t p = 0)
        : type(t), pos(p) {}

    static
    char
    get_type_label(const map_t t) {
        switch (t) {
        case NONE: return 'N';
        case MATCH: return 'M';
        case INSERT: return 'I';
        case SOFT_CLIP: return 'S';
        case CONFLICT: return 'X';
        default:
            assert(0);
            return '\0';
        }
    }

    map_t type;
    pos_t pos;
};

}


#ifdef DEBUG_ALIGN_CLIP
static
void
dump_ref_map(const std::vector<ref_map_type>& ref_map,
             std::ostream& os) {

    const unsigned rs(ref_map.size());
    for (unsigned i(0); i<rs; ++i) {
        const ref_map_type& rmt(ref_map[i]);
        os << "rm: " << i << " " << ref_map_type::get_type_label(rmt.type) << " " << rmt.pos << "\n";
    }
}
#endif



static
void
get_alignment_ref_map(const alignment& al,
                      std::vector<ref_map_type>& ref_map) {

    using namespace ALIGNPATH;

    ref_map.clear();

    pos_t ref_head_pos(al.pos);

    const unsigned as(al.path.size());
    for (unsigned i(0); i<as; ++i) {
        const path_segment& ps(al.path[i]);

        if       (ps.type == MATCH) {
            for (unsigned j(0); j<ps.length; ++j) {
                ref_map.push_back(ref_map_type(ref_map_type::MATCH,ref_head_pos+j));
            }
            ref_head_pos += ps.length;
        } else if (ps.type == INSERT) {
            for (unsigned j(0); j<ps.length; ++j) {
                ref_map.push_back(ref_map_type(ref_map_type::INSERT,0));
            }
        } else if ((ps.type == DELETE) or (ps.type == SKIP)) {
            ref_head_pos += ps.length;
        } else if (ps.type == SOFT_CLIP) {
            for (unsigned j(0); j<ps.length; ++j) {
                ref_map.push_back(ref_map_type(ref_map_type::SOFT_CLIP,0));
            }
        } else if (ps.type == HARD_CLIP) {
            // do nothing...
        } else {
            std::ostringstream oss;
            oss << "Can't handle cigar code: " << segment_type_to_cigar_code(ps.type) << "\n";
            throw blt_exception(oss.str().c_str());
        }
    }
}



static
void
mark_ref_map_conflicts(const alignment& al,
                       std::vector<ref_map_type>& ref_map) {

    using namespace ALIGNPATH;

    pos_t ref_head_pos(al.pos);
    pos_t read_head_pos(0);

    const unsigned as(al.path.size());
    for (unsigned i(0); i<as; ++i) {
        const path_segment& ps(al.path[i]);

        if       (ps.type == MATCH) {
            for (unsigned j(0); j<ps.length; ++j) {
                ref_map_type& rm(ref_map[read_head_pos+j]);
                if (rm.type == ref_map_type::CONFLICT) continue;
                if ((rm.type != ref_map_type::MATCH) or
                    (rm.pos != (ref_head_pos+static_cast<pos_t>(j)))) {
                    rm.type=ref_map_type::CONFLICT;
                }
            }
            read_head_pos += ps.length;
            ref_head_pos += ps.length;
        } else if (ps.type == INSERT) {
            for (unsigned j(0); j<ps.length; ++j) {
                ref_map_type& rm(ref_map[read_head_pos+j]);
                if (rm.type == ref_map_type::CONFLICT) continue;
                if (rm.type != ref_map_type::INSERT) {
                    rm.type=ref_map_type::CONFLICT;
                }
            }
            read_head_pos += ps.length;
        } else if ((ps.type == DELETE) or (ps.type == SKIP)) {
            ref_head_pos += ps.length;
        } else if (ps.type == SOFT_CLIP) {
            for (unsigned j(0); j<ps.length; ++j) {
                ref_map_type& rm(ref_map[read_head_pos+j]);
                if (rm.type == ref_map_type::CONFLICT) continue;
                if (rm.type != ref_map_type::SOFT_CLIP) {
                    rm.type=ref_map_type::CONFLICT;
                }
            }
            read_head_pos += ps.length;
        } else if (ps.type == HARD_CLIP) {
            // do nothing...
        } else {
            std::ostringstream oss;
            oss << "Can't handle cigar code: " << segment_type_to_cigar_code(ps.type) << "\n";
            throw blt_exception(oss.str().c_str());
        }
    }
}



static
void
extend_or_add_sc(alignment& al,
                 const unsigned length) {
    using namespace ALIGNPATH;

    if ((! al.path.empty()) && (al.path.back().type == SOFT_CLIP)) {
        al.path.back().length += length;
    } else {
        al.path.push_back(path_segment(SOFT_CLIP,length));
    }
}



static
void
soft_clip_alignment(alignment& al,
                    const unsigned leading_clip,
                    const unsigned trailing_clip) {

    using namespace ALIGNPATH;

    assert(leading_clip < trailing_clip);

    unsigned read_head_pos(0);

    alignment new_al;
    new_al.pos=al.pos;
    new_al.is_fwd_strand=al.is_fwd_strand;

    const unsigned as(al.path.size());
    for (unsigned i(0); i<as; ++i) {
        path_segment& ps(al.path[i]);

        if ((ps.type == MATCH) or (ps.type == INSERT)) {
            if (leading_clip > read_head_pos) {
                const unsigned clip_length(std::min(ps.length,(leading_clip-read_head_pos)));
                extend_or_add_sc(new_al,clip_length);
                if (ps.type == MATCH) new_al.pos += static_cast<pos_t>(clip_length);
                if (clip_length < ps.length) {
                    const unsigned frac_length(ps.length-clip_length);
                    new_al.path.push_back(path_segment(ps.type,frac_length));
                }
            } else if (trailing_clip < (read_head_pos+ps.length)) {
                const unsigned clip_length(std::min(ps.length,((read_head_pos+ps.length)-trailing_clip)));
                if (clip_length < ps.length) {
                    const unsigned frac_length(ps.length-clip_length);
                    new_al.path.push_back(path_segment(ps.type,frac_length));
                }
                extend_or_add_sc(new_al,clip_length);
            } else {
                new_al.path.push_back(ps);
            }
            read_head_pos += ps.length;
        } else if ((ps.type == DELETE) or (ps.type == SKIP)) {
            if (leading_clip >= read_head_pos) {
                new_al.pos += static_cast<pos_t>(ps.length);
            } else if (trailing_clip <= read_head_pos) {
            } else {
                new_al.path.push_back(ps);
            }
        } else if (ps.type == SOFT_CLIP) {
            extend_or_add_sc(new_al,ps.length);
            read_head_pos += ps.length;
        } else if (ps.type == HARD_CLIP) {
            new_al.path.push_back(ps);
        } else {
            std::ostringstream oss;
            oss << "Can't handle cigar code: " << segment_type_to_cigar_code(ps.type) << "\n";
            throw blt_exception(oss.str().c_str());
        }
    }

    al=new_al;
}




void
get_clipped_alignment_from_cal_pool(const cal_pool_t& max_cal_pool,
                                    const unsigned best_cal_id,
                                    alignment& al) {

    const unsigned n_cal(max_cal_pool.size());
    assert(n_cal>0);
    assert(best_cal_id<n_cal);

    al=max_cal_pool[best_cal_id]->al;
    if (n_cal==1) return;

    // create read_pos->ref_pos map for first alignment -- start
    // marking off the conflict positions in other alignments:
    std::vector<ref_map_type> ref_map;
    get_alignment_ref_map(al,ref_map);
#ifdef DEBUG_ALIGN_CLIP
    std::cerr << "VARMIT: initial ref_map:\n";
    dump_ref_map(ref_map,std::cerr);
#endif
    for (unsigned i(0); i<n_cal; ++i) {
        if (i==best_cal_id) continue;
        mark_ref_map_conflicts(max_cal_pool[i]->al,ref_map);
#ifdef DEBUG_ALIGN_CLIP
        std::cerr << "VARMIT: modified ref_map round " << i << ":\n";
        dump_ref_map(ref_map,std::cerr);
#endif
    }

    // an internal ambiguous alignment path might be possible, but
    // very rare, so we'll ignore anything which doesn't occur on
    // the ends:
    //
    // from each end -- move in until a match is found, then move out
    // until a conflict (or the end) is found:
    const unsigned read_size(ref_map.size());

    unsigned leading_clip(0);
    for (; leading_clip<read_size; leading_clip++) {
        if (ref_map[leading_clip].type == ref_map_type::MATCH) break;
    }
    for (; leading_clip>0; leading_clip--) {
        if ((ref_map[leading_clip-1].type == ref_map_type::CONFLICT) ||
            (ref_map[leading_clip-1].type == ref_map_type::SOFT_CLIP)) break;
    }

    unsigned trailing_clip(read_size);
    for (; trailing_clip>0; trailing_clip--) {
        if (ref_map[trailing_clip-1].type == ref_map_type::MATCH) break;
    }
    for (; trailing_clip<read_size; trailing_clip++) {
        if ((ref_map[trailing_clip].type == ref_map_type::CONFLICT) ||
            (ref_map[trailing_clip].type == ref_map_type::SOFT_CLIP)) break;
    }
#ifdef DEBUG_ALIGN_CLIP
    std::cerr << "VARMIT: lead,trail: " << leading_clip << " " << trailing_clip << "\n";
#endif

    if (leading_clip>=trailing_clip) {
        al.clear();
        return;
    }

    if ((leading_clip!=0) ||
        (trailing_clip!=read_size)) {
        soft_clip_alignment(al,leading_clip,trailing_clip);
    }
}

