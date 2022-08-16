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
#ifndef __NPLOID_GENOTYPE_UTIL_HH
#define __NPLOID_GENOTYPE_UTIL_HH

#include "blt_util/seq_util.hh"

#include <boost/utility.hpp>

#include <string>
#include <vector>


unsigned
nploid_gtype_size(const unsigned ploidy);


/// pre-compute genotype count, labels and allele frequencies for any
/// given ploidy:
///
struct nploid_info : private boost::noncopyable {

    nploid_info(const unsigned init_ploidy);

    unsigned
    ploidy() const { return _ploidy; }

    unsigned
    gtype_size() const { return _ginfo.size(); }

    const char*
    label(const unsigned gt_index) const {
        return _ginfo[gt_index].label.c_str();
    }

    unsigned
    get_ref_gtype(const char a) const {
        return _ref_gtype[base_to_id(a)];
    }

    /// look up expectation frequencies directly:
    double
    expect_freq(const unsigned gt_index,
                const unsigned base_id) const {
        return expect_freq_level(gt_index,base_id)*expect_freq_chunk();
    }

    /// expectation frequencies will always be multiples of 1/ploidy,
    /// so more efficient client code can be written by asking for freq
    /// level and chunk size, where expect=level*chunk, chunk=1/ploidy
    ///
    double
    expect_freq_chunk() const { return _echunk; }

    unsigned
    expect_freq_level(const unsigned gt_index,
                      const unsigned base_id) const {
        return _ginfo[gt_index].efreq_levels[base_id];
    }

    unsigned
    expect_freq_level_size() const { return _ploidy+1; }

private:

    struct gtype_info {
        gtype_info() {
            for (unsigned i(0); i<N_BASE; ++i) efreq_levels[i]=0;
        }

        gtype_info(const gtype_info& g) : label(g.label) {
            for (unsigned i(0); i<N_BASE; ++i) efreq_levels[i]=g.efreq_levels[i];
        }

        gtype_info& operator=(const gtype_info& rhs) {
            if ( this == &rhs ) return *this;
            label=rhs.label;
            for (unsigned i(0); i<N_BASE; ++i) efreq_levels[i]=rhs.efreq_levels[i];
            return *this;
        }

        std::string label;
        unsigned efreq_levels[N_BASE];
    };

    void
    ginfo_init(const unsigned init_ploidy,
               const unsigned pli,
               const unsigned init_i,
               gtype_info& gi);

    void
    ref_gtype_init(const unsigned init_ploidy,
                   const unsigned pli,
                   const unsigned init_i,
                   const bool is_hom,
                   unsigned*& ref,
                   unsigned& index);

    unsigned _ploidy;
    double _echunk;
    std::vector<gtype_info> _ginfo;
    unsigned _ref_gtype[N_BASE];
};


#endif
