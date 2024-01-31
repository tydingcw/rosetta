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
#include <core/kinematics/FoldTree.hh>

/// Project headers
//#include <core/types.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/moves/DsspMover.hh>

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

    core::kinematics::FoldTree fold_tree_from_dssp_string(
            std::string input_string
            //utility::vector1< std::pair< core::Size, core::Size > > &input_vect
    ) {
        auto input_vect = identify_secondary_structure_spans(input_string);
        std::cout << "New string length " << input_string.size() << std::endl;
        std::cout << "New vector length " << input_vect.size() << std::endl;

        core::kinematics::FoldTree ft;
        //core::Size edges = 4*(input_vect.size()-2);
        //core::Size jumps = 2*(input_vect.size()-2);

        //Making the start
        std::pair< core::Size, core::Size > node1 = input_vect[1];
        core::Size center = (node1.first + node1.second)/2;
        std::cout << "New edge " << center << " " << 1 << std::endl;
        std::cout << "New edge " << center << " " << node1.second << std::endl;
        ft.add_edge( center, 1, core::kinematics::Edge::PEPTIDE );
        ft.add_edge( center, node1.second, core::kinematics::Edge::PEPTIDE );

        core::Size count = 1;
        //Creating the foldtree for given SS and previous loop
        for (size_t i = 2; i <= input_vect.size() - 1; i++){
            //assuming loop size 3+ between SS elements
            std::cout << "starting " << i << std::endl;
            std::pair< core::Size, core::Size > node_i = input_vect[i];
            core::Size center_i = (node_i.first + node_i.second)/2;

            std::pair< core::Size, core::Size > node_prev = input_vect[i-1];
            core::Size loop_start = node_prev.second + 1; //prev node is SS, need loop
            core::Size loop_end = node_i.first -1;
            core::Size loop_center = (loop_start + loop_end)/2;

            std::cout << "New jump " << center << " " << loop_center << " " << count << std::endl;
            std::cout << "New edge " << loop_center << " " << loop_start << std::endl;
            std::cout << "New edge " << loop_center << " " << loop_end << std::endl;
            ft.add_edge( center, loop_center, count );
            ft.add_edge( loop_center, loop_start, core::kinematics::Edge::PEPTIDE );
            ft.add_edge( loop_center, loop_end, core::kinematics::Edge::PEPTIDE );

            count++;

            std::cout << "New jump " << center << " " << center_i << " " << count << std::endl;
            std::cout << "New edge " << center_i << " " << node_i.first << std::endl;
            std::cout << "New edge " << center_i << " " << node_i.second << std::endl;
            ft.add_edge( center, center_i, count );
            ft.add_edge( center_i, node_i.first, core::kinematics::Edge::PEPTIDE );
            ft.add_edge( center_i, node_i.second, core::kinematics::Edge::PEPTIDE );

            count++;
        }

        //Making the end
        std::pair< core::Size, core::Size > node_end = input_vect[input_vect.size()];

        std::pair< core::Size, core::Size > node_prev = input_vect[input_vect.size()-1];
        core::Size loop_start = node_prev.second + 1; //prev node is SS, need loop
        core::Size loop_end = node_end.first -1;
        core::Size loop_center = (loop_start + loop_end)/2;

        std::cout << "New jump " << center << " " << loop_center << " " << count << std::endl;
        std::cout << "New edge " << loop_center << " " << loop_start << std::endl;
        std::cout << "New edge " << loop_center << " " << loop_end << std::endl;
        ft.add_edge( center, loop_center, count );
        ft.add_edge( loop_center, loop_start, core::kinematics::Edge::PEPTIDE );
        ft.add_edge( loop_center, loop_end, core::kinematics::Edge::PEPTIDE );
        count++;

        core::Size center_end = (node_end.first + node_end.second)/2;
        std::cout << "New jump " << center << " " << center_end << " " << count << std::endl;
        std::cout << "New edge " << center_end << " " << node_end.first << std::endl;
        std::cout << "New edge " << center_end << " " << input_string.size() << std::endl;
        ft.add_edge( center, center_end, count );
        ft.add_edge( center_end, node_end.first, core::kinematics::Edge::PEPTIDE );
        //Need to get last element
        ft.add_edge( center_end, input_string.size(), core::kinematics::Edge::PEPTIDE );

        return ft;
    }

    core::kinematics::FoldTree fold_tree_from_ss(core::pose::Pose & pose) {
        core::scoring::dssp::Dssp dssp_struct = core::scoring::dssp::Dssp(pose);
        //core::scoring::dssp::Dssp dssp_mem = core::scoring::dssp::Dssp();
        std::string input_string = dssp_struct.get_dssp_secstruct();
        return fold_tree_from_dssp_string(input_string);
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

    void test_fold_tree_size() {
        auto my_tree = fold_tree_from_dssp_string("   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     ");
        std::cout << "Size " << my_tree.size() << std::endl;
        TS_ASSERT( my_tree.size() == 38 );
		}

    void test_pdb_fold_tree() {
        core::pose::Pose test_pose = create_test_in_pdb_pose();
        auto my_tree = fold_tree_from_ss( test_pose );
        //std::cout << "Size " << my_tree.size() << std::endl;
        test_pose.fold_tree(my_tree);
        TS_ASSERT( my_tree.check_fold_tree() );
    }

    void test_empty() {
        TS_ASSERT( identify_secondary_structure_spans( "" ).size() == 0 );
    }

    void test1_identify_secondary_structure_spans() {
        //7 secondary structure elements, spanning residues
        // 4 to 8, 12 to 19, 22 to 26, 36 to 41, 45 to 55, 58 to 62, and 65 to 68
        std::string test_1 = "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   ";
        auto output_1 = identify_secondary_structure_spans( test_1 );
        utility::vector1< std::pair< core::Size, core::Size > > result_1 = {
                {4,8},
                {12,19},
                {22, 26},
                {36, 41},
                {45,55},
                {58,62},
                {65, 68}
        };
        TS_ASSERT( output_1 == result_1 );
    }

    void test2_identify_secondary_structure_spans() {
        //7 secondary structure elements, spanning residues
        // 4 to 8, 12 to 19, 22 to 26, 36 to 41, 45 to 55, 58 to 62, and 65 to 68
        std::string test_1 = "HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH ";
        auto output_1 = identify_secondary_structure_spans( test_1 );
        utility::vector1< std::pair< core::Size, core::Size > > result_1 = {
                {1,  7},
                {11, 22},
                {29, 40},
                {41, 50},
                {51, 57},
                {59, 62},
                {63, 65}
        };
        TS_ASSERT( output_1 == result_1 );
    }

    void test3_identify_secondary_structure_spans() {
        //7 secondary structure elements, spanning residues
        // 4 to 8, 12 to 19, 22 to 26, 36 to 41, 45 to 55, 58 to 62, and 65 to 68
        std::string test_1 = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE";
        auto output_1 = identify_secondary_structure_spans( test_1 );
        utility::vector1< std::pair< core::Size, core::Size > > result_1 = {
                {1,  9},
                {11, 18},
                {20, 28},
                {30, 30},
                {32, 36},
                {38, 38},
                {40, 40},
                {42, 42},
                {44, 51}
        };
        TS_ASSERT( output_1 == result_1 );
    }

};
