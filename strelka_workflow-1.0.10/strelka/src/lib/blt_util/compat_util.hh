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
// take care of some (mostly C99) functions not available in VS C++
//

#ifndef __COMPAT_UTIL
#define __COMPAT_UTIL


#include <string>


#ifdef _WIN32
#define snprintf _snprintf
#endif


double
compat_round(const double x);


const char*
compat_basename(const char* s);


// gets canonical name of paths, but only when these refer to existing items
// returns false on error.
bool
compat_realpath(std::string& path);


#endif
