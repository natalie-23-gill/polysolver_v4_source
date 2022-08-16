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


#ifndef __CHROM_DEPTH_MAP
#define __CHROM_DEPTH_MAP


#include <map>
#include <string>


typedef std::map<std::string,double> cdmap_t;


// parse the chrom depth file
void
parse_chrom_depth(const std::string& chrom_depth_file,
                  cdmap_t& chrom_depth);


#endif
