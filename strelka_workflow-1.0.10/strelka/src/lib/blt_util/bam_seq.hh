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

/// \author Chris Saunders
///
#ifndef __BAM_SEQ_HH
#define __BAM_SEQ_HH

#include "blt_util/blt_types.hh"
#include "blt_util/reference_contig_segment.hh"

#include <stdint.h>

#include <cassert>
#include <ciso646>
#include <iosfwd>
#include <string>


namespace BAM_BASE {
enum index_t {
    REF = 0x0,
    A = 0x1,
    C = 0x2,
    G = 0x4,
    T = 0x8,
    ANY = 0xF
};
}


inline
char
get_bam_seq_char(const uint8_t a) {

    using namespace BAM_BASE;

    switch (a) {
    case REF: return '=';
    case A:   return 'A';
    case C:   return 'C';
    case G:   return 'G';
    case T:   return 'T';
    default:  return 'N';
    }
}


inline
uint8_t
get_bam_seq_code(const char c) {

    using namespace BAM_BASE;

    switch (c) {
    case '=': return REF;
    case 'A': return A;
    case 'C': return C;
    case 'G': return G;
    case 'T': return T;
    default:  return ANY;
    }
}


// interface to bam_seq -- allows us to pass either compressed
// sequences from bam files and regular strings using the same
// object:
//
struct bam_seq_base {

    virtual ~bam_seq_base() {}

    virtual uint8_t get_code(pos_t i) const = 0;

    virtual char get_char(const pos_t i) const = 0;

    virtual unsigned size() const = 0;

protected:
    bool
    is_in_range(const pos_t i) const {
        return ((i>=0) && (i<static_cast<pos_t>(size())));
    }
};

std::ostream& operator<<(std::ostream& os, const bam_seq_base& bs);


//
//
struct bam_seq : public bam_seq_base {

    bam_seq(const uint8_t* s,
            const uint16_t init_size,
            const uint16_t offset=0)
        : _s(s), _size(init_size), _offset(offset) {}

#if 0
    bam_seq(const bam_seq bs,
            const uint16_t size,
            const uint16_t offset=0)
        : _s(bs.s), _size(size), _offset(bs.offset+offset) {
        assert((offset+size)<=bs.size);
    }
#endif

    uint8_t
    get_code(pos_t i) const {
        if (! is_in_range(i)) return BAM_BASE::ANY;
        i += static_cast<pos_t>(_offset);
        return _s[(i/2)] >> 4*(1-(i%2)) & 0xf;
    }

    char
    get_char(const pos_t i) const {
        return get_bam_seq_char(get_code(i));
    }

    unsigned size() const { return _size; }

private:
    const uint8_t* _s;
    uint16_t _size;
    uint16_t _offset;
};


//
//
struct string_bam_seq : public bam_seq_base {

    string_bam_seq(const std::string& s)
        : _s(s.c_str()), _size(s.size()) {}

    string_bam_seq(const char* s,
                   const unsigned init_size)
        : _s(s), _size(init_size) {}

    uint8_t
    get_code(pos_t i) const {
        return get_bam_seq_code(get_char(i));
    }

    char
    get_char(const pos_t i) const {
        if (! is_in_range(i)) return 'N';
        return _s[i];
    }

    unsigned size() const { return _size; }

private:
    const char* _s;
    unsigned _size;
};


//
//
struct rc_segment_bam_seq : public bam_seq_base {

    rc_segment_bam_seq(const reference_contig_segment& r)
        : _r(r)
    {}

    uint8_t
    get_code(pos_t i) const {
        return get_bam_seq_code(get_char(i));
    }

    char
    get_char(const pos_t i) const {
        return _r.get_base(i);
    }

    unsigned size() const { return _r.end(); }

private:
    const reference_contig_segment& _r;
};


#endif
