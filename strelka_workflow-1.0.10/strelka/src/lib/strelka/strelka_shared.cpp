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

#include "position_somatic_snv.hh"
#include "position_somatic_snv_grid.hh"
#include "position_somatic_snv_strand_grid.hh"
#include "somatic_indel.hh"
#include "somatic_indel_grid.hh"

#include "strelka/strelka_shared.hh"



strelka_deriv_options::
strelka_deriv_options(const strelka_options& opt,
                      const reference_contig_segment& ref)
    : base_t(opt,ref)
    , _sscaller_strand_grid(new somatic_snv_caller_strand_grid(opt,pdcaller()))
    , _sicaller_grid(new somatic_indel_caller_grid(opt,incaller()))
{}



strelka_deriv_options::
~strelka_deriv_options() {}
