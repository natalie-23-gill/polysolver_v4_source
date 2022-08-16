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
#ifndef __SIG_HANDLER_HH
#define __SIG_HANDLER_HH

/// \brief logs sigint/sigterm events (with cmdline):
///
void
initialize_blt_signals(const char* progname,
                       const char* cmdline);

#endif
