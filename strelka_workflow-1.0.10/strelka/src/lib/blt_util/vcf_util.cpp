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

#include "blt_exception.hh"
#include "vcf_util.hh"

#include <cassert>
#include <cctype>
#include <ctime>

#include <iostream>
#include <sstream>



std::ostream&
vcf_fileDate(std::ostream& os) {
    const time_t t(time(NULL));
    struct tm* ct(localtime(&t));
    assert(NULL != ct);

    static const unsigned dsize(64);
    char datebuf[dsize];
    const size_t ret(strftime(datebuf,dsize,"%Y%m%d",ct));
    assert(ret!=0);
    return os << datebuf;
}



struct gt_parse_helper {

    // return is_valid_genotype
    static
    bool
    start(const char* gt,
          std::vector<int>& gti,
          const bool is_badend) {
        gti.clear();
        if (isdigit(*gt)) return digit(gt,gti,is_badend);

        switch (*gt) {
        case '.' :  return unknown(gt,gti,is_badend);
        default: return false;
        }
    }

private:

    static
    bool
    unknown(const char* gt,
            std::vector<int>& gti,
            const bool is_badend) {
        gt++;
        gti.push_back(-1);
        switch (*gt) {
        case '\0' : return true;
        case '|' :
        case '/' : return sep(gt,gti,is_badend);
        default : return is_badend;
        }
    }

    static
    bool
    sep(const char* gt,
        std::vector<int>& gti,
        const bool is_badend) {
        gt++;
        if (isdigit(*gt)) return digit(gt,gti,is_badend);
        switch (*gt) {
        case '.' : return unknown(gt,gti,is_badend);
        default : return false;
        }
    }

    static
    bool
    digit(const char* gt,
          std::vector<int>& gti,
          const bool is_badend) {
        int val(0);
        while (isdigit(*gt)) {
            val = val*10 + static_cast<int>(*gt-'0');
            gt++;
        }
        gti.push_back(val);

        switch (*gt) {
        case '\0' : return true;
        case '|' :
        case '/' : return sep(gt,gti,is_badend);
        default : return is_badend;
        }
    }
};




void
parse_gt(const char* gt,
         std::vector<int>& gti,
         const bool is_allow_bad_end_char) {

    if (! gt_parse_helper::start(gt,gti,is_allow_bad_end_char)) {
        std::ostringstream oss;
        oss << "ERROR: can't parse genotype string: '" << gt << "'\n";
        throw blt_exception(oss.str().c_str());
    }
}
