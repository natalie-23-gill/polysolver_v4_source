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

#ifndef __STRELKA_INFO_HH
#define __STRELKA_INFO_HH

#include "blt_util/prog_info.hh"


struct strelka_info : public prog_info {

    static
    const prog_info& get() {
        static const strelka_info vci;
        return vci;
    }

private:
    const char* name() const {
        static const char NAME[] = "strelka";
        return NAME;
    }

    const char* version() const {
        static const char VERSION[] = "2.0.12";
        return VERSION;
    }

    void usage(const char* xmessage = 0) const;

    void doc() const;

    strelka_info() {}
    virtual ~strelka_info() {}
};

#endif
