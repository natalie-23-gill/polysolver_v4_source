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

#ifndef __STARLING_ARG_PARSE_H
#define __STARLING_ARG_PARSE_H

#include "blt_common/blt_arg_parse_util.hh"
#include "starling_common/starling_shared.hh"


// this combines starling argument parsing and validation:
//
void
legacy_starling_arg_parse(arg_data& ad,
                          starling_options& client_opt);

#endif
