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

#include "blt_util/bam_util.hh"
#include "blt_util/bam_seq.hh"
#include "blt_util/log.hh"

#include <cassert>

#include <iostream>



static
void
change_bam_data_len(const int new_len,
                    bam1_t& br) {

    assert(new_len>=0);

    if (new_len > br.m_data) {
        br.m_data = new_len;
        kroundup32(br.m_data);
        br.data = (uint8_t*) realloc(br.data,br.m_data);
        if (NULL == br.data) {
            log_os << "ERROR: failed to realloc BAM data size to: " << new_len << "\n";
            exit(EXIT_FAILURE);
        }
    }
    br.data_len = new_len;
}



void
change_bam_data_segment_len(const int end,
                            const int delta,
                            bam1_t& br) {

    assert(end>=0);
    if (0==delta) return;
    const int old_len(br.data_len);
    const int new_len(old_len+delta);
    const int tail_size(old_len-end);
    assert(tail_size>=0);
    change_bam_data_len(new_len,br);

    // move post-segment data to its new position:
    if (0==tail_size) return;
    uint8_t* old_tail_ptr(br.data+end);
    uint8_t* new_tail_ptr(old_tail_ptr+delta);
    memmove(new_tail_ptr,old_tail_ptr,tail_size);
}



void
edit_bam_qname(const char* name,
               bam1_t& br) {

    bam1_core_t& bc(br.core);

    const uint32_t tmp_size(strlen(name)+1);
    if (tmp_size & 0xffffff00) {
        log_os << "ERROR: name is too long to be entered in BAM qname field: " << name << "\n";
        exit(EXIT_FAILURE);
    }
    const uint8_t new_qname_size(tmp_size);
    const uint8_t old_qname_size(bc.l_qname);
    const int delta(new_qname_size-old_qname_size);

    if (0 != delta) {
        change_bam_data_segment_len(old_qname_size,delta,br);
        bc.l_qname=new_qname_size;
    }

    strcpy(bam1_qname(&br),name);
}


static
inline
int seq_size(const int a) { return a+(a+1)/2; }


void
edit_bam_read_and_quality(const char* read,
                          const uint8_t* qual,
                          bam1_t& br) {

    const int new_len(strlen(read));
    const int old_size(seq_size(br.core.l_qseq));
    const int new_size(seq_size(new_len));
    const int delta(new_size-old_size);

    if (0 != delta) {
        const int end(bam1_aux(&br)-br.data);
        change_bam_data_segment_len(end,delta,br);
    }
    br.core.l_qseq = new_len;
    // update seq:
    uint8_t* p(bam1_seq(&br));
    memset(p,0,(new_len+1)/2);
    for (int i(0); i<new_len; ++i) {
        p[i/2] |= get_bam_seq_code(read[i]) << 4*(1-i%2);
    }
    // update qual
    memcpy(bam1_qual(&br),qual,new_len);
}



void
nuke_bam_aux_field(bam1_t& br,
                   const char* tag) {

    while (true) {
        uint8_t* p(bam_aux_get(&br,tag));
        if (NULL==p) return;
        bam_aux_del(&br,p);
    }
}



// optimize storage of an unsigned int in bam
void
bam_aux_append_unsigned(bam1_t& br,
                        const char* tag,
                        uint32_t x) {

    if       (x & 0xffff0000) {
        bam_aux_append(&br,tag,'I',4,reinterpret_cast<uint8_t*>(&x));
    } else if (x & 0xff00) {
        uint16_t y(x);
        bam_aux_append(&br,tag,'S',2,reinterpret_cast<uint8_t*>(&y));
    } else {
        uint8_t z(x);
        bam_aux_append(&br,tag,'C',1,&z);
    }
}



bool
check_header_compatibility(const bam_header_t* h1,
                           const bam_header_t* h2) {

    if (h1->n_targets != h2->n_targets) {
        return false;
    }

    for (int32_t i(0); i<h1->n_targets; ++i) {
        if (h1->target_len[i] != h2->target_len[i]) return false;
        if (0 != strcmp(h1->target_name[i],h2->target_name[i])) return false;
    }
    return true;
}

