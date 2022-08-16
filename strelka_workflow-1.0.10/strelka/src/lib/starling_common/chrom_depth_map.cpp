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


#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/parse_util.hh"
#include "starling_common/chrom_depth_map.hh"

#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>



// parse the chrom depth file
void
parse_chrom_depth(const std::string& chrom_depth_file,
                  cdmap_t& chrom_depth) {

    if (chrom_depth_file.empty()) return;

    std::ifstream depth_is(chrom_depth_file.c_str());
    if (! depth_is) {
        log_os << "ERROR: Failed to open chrom depth file '" << chrom_depth_file << "'\n";
        exit(EXIT_FAILURE);
    }

    static const unsigned buff_size(1024);
    char buff[buff_size];

    unsigned line_no(0);

    while (true) {
        depth_is.getline(buff,buff_size);
        if (! depth_is) {
            if     (depth_is.eof()) break;
            else {
                log_os << "ERROR: unexpected failure while attempting to read chrom depth file line " << (line_no+1) << "\n";
                exit(EXIT_FAILURE);
            }
        } else {
            ++line_no;
        }

        char* word2(strchr(buff,'\t'));
        if (NULL == word2) {
            log_os << "ERROR: unexpected format in read chrom depth file line " << (line_no) << "\n";
            exit(EXIT_FAILURE);
        }
        *(word2++) = '\0';
        try {
            const char* s(word2);
            chrom_depth[buff] = illumina::blt_util::parse_double(s);
        } catch (const blt_exception& e) {
            log_os << "ERROR: unexpected format in read chrom depth file line " << (line_no) << "\n";
            throw;
        }
    }
}

