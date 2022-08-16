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

#ifndef __ALIGN_PATH_BAM_UTIL_HH
#define __ALIGN_PATH_BAM_UTIL_HH

#include "blt_util/bam_util.hh"
#include "starling_common/align_path.hh"


// convert internal BAM cigar representation directly into a path:
//
void
bam_cigar_to_apath(const uint32_t* bam_cigar,
                   const unsigned n_cigar,
                   ALIGNPATH::path_t& apath);

// convert apath to internal BAM cigar representation:
//
// bam_cigar should already be set to apath.size() capacity
//
void
apath_to_bam_cigar(const ALIGNPATH::path_t& apath,
                   uint32_t* bam_cigar);

// convert apath into a CIGAR string and replace CIGAR in BAM record
//
void
edit_bam_cigar(const ALIGNPATH::path_t& apath,
               bam1_t& br);

#endif
