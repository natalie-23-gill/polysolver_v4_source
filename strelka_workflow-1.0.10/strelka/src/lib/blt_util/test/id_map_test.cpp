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

#include "boost/test/unit_test.hpp"

#include "id_map.hh"

#include <string>


BOOST_AUTO_TEST_SUITE( test_id_map )


BOOST_AUTO_TEST_CASE( test_id_set ) {

    id_set<std::string> iset;

    iset.insert_key("brown");
    iset.insert_key("fox");

    BOOST_CHECK_EQUAL(iset.test_key("brown"),true);
    BOOST_CHECK_EQUAL(iset.test_key("x"),false);

    const unsigned expect_id(1);

    BOOST_CHECK_EQUAL(iset.get_id("fox"),expect_id);
    BOOST_CHECK_EQUAL(iset.get_key(expect_id),std::string("fox"));
}


BOOST_AUTO_TEST_CASE( test_id_map ) {

    id_map<std::string,std::string> imap;

    imap.insert("brown","123");
    imap.insert("fox","456");

    BOOST_CHECK_EQUAL(imap.test_key("brown"),true);
    BOOST_CHECK_EQUAL(imap.test_key("x"),false);

    const unsigned expect_id(1);

    BOOST_CHECK_EQUAL(imap.get_id("fox"),expect_id);
    BOOST_CHECK_EQUAL(imap.get_key(expect_id),std::string("fox"));
    BOOST_CHECK_EQUAL(imap.get_value(expect_id),std::string("456"));
}

BOOST_AUTO_TEST_SUITE_END()

