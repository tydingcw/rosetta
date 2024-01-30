// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

//#include <protocols/match/upstream/ProteinSCSampler.hh>
//#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers

/// Project headers
//#include <core/types.hh>

// C++ headers
#include <utility/vector1.hh>
#include <string>
#include <iostream>

//Auto Headers
//#include <core/pack/dunbrack/DunbrackRotamer.hh>


//using namespace protocols::match;
//using namespace protocols::match::upstream;


// --------------- Test Class --------------- //

class FoldTreeFromSSTests : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

    utility::vector1< std::pair< core::Size, core::Size > >
    identify_secondary_structure_spans( std::string const & ss_string )
    {
        utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
        core::Size strand_start = -1;
        for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
            if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
                if ( int( strand_start ) == -1 ) {
                    strand_start = ii;
                } else if ( ss_string[ii] != ss_string[strand_start] ) {
                    ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
                    strand_start = ii;
                }
            } else {
                if ( int( strand_start ) != -1 ) {
                    ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
                    strand_start = -1;
                }
            }
        }
        if ( int( strand_start ) != -1 ) {
            // last residue was part of a ss-eleemnt
            ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
        }
        for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
            std::cout << "SS Element " << ii << " from residue "
                      << ss_boundaries[ ii ].first << " to "
                      << ss_boundaries[ ii ].second << std::endl;
        }
        return ss_boundaries;
    }


    // Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	void test_hello_world() {
		TS_ASSERT( true );
	}

    void test_empty() {
        TS_ASSERT( identify_secondary_structure_spans( "" ).size() == 0 );
    }

    //void test_identify_secondary_structure_spans() {
        //7 secondary structure elements, spanning residues
        // 4 to 8, 12 to 19, 22 to 26, 36 to 41, 45 to 55, 58 to 62, and 65 to 68
    //    std::string test_1 = "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   ";
    //}
};
