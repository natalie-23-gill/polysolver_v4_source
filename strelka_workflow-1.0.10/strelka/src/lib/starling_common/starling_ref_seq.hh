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

#ifndef __STARLING_REF_SEQ_HH
#define __STARLING_REF_SEQ_HH

#include "starling_common/starling_shared.hh"

#include <string>


void
get_starling_ref_seq(const starling_options& opt,
                     reference_contig_segment& ref);


#endif
