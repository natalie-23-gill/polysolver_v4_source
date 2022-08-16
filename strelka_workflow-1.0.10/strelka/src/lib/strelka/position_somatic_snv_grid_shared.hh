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

/// variation on the original strowman snv caller -- implements a
/// compile-time specified grid in allele frequency space and requires
/// similar frequency as definition of non-somatic.
///

/// \author Chris Saunders
///
#ifndef __POSITION_SOMATIC_SNV_GRID_SHARED_HH
#define __POSITION_SOMATIC_SNV_GRID_SHARED_HH


struct somatic_snv_genotype_grid {

    somatic_snv_genotype_grid()
        : is_snv(false) {}

    typedef bool tier_t;

    struct result_set {
        unsigned ntype;
        unsigned max_gt;
        int snv_qphred;
        int snv_from_ntype_qphred;
    };

    bool is_snv;
    tier_t snv_tier;
    tier_t snv_from_ntype_tier;
    unsigned ref_gt;
    result_set rs;
};

#endif
