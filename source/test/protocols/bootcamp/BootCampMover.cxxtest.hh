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
#include <protocols/moves/MoverFactory.hh>
#include <basic/datacache/DataMap.hh>
//#include <basic/resource_manager/ResourceManager.hh>
#include <test/util/rosettascripts.hh>

/// Project headers
//#include <core/types.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/bootcamp/BootCampMover.hh>

// C++ headers
#include <utility/vector1.hh>
#include <string>
#include <iostream>

//Auto Headers
//#include <core/pack/dunbrack/DunbrackRotamer.hh>


//using namespace protocols::match;
//using namespace protocols::match::upstream;
using namespace protocols::bootcamp;


// --------------- Test Class --------------- //

class BootCampMoverTest : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.


    // Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

//    void test_pdb_fold_tree() {
//        core::pose::Pose test_pose = create_test_in_pdb_pose();
//        auto my_tree = fold_tree_from_ss( test_pose );
//        //std::cout << "Size " << my_tree.size() << std::endl;
//        test_pose.fold_tree(my_tree);
//        TS_ASSERT( my_tree.check_fold_tree() );
//    }

    void test_set_num_iterations() {
        core::Size iters = 100;
        auto boot_mv = protocols::bootcamp::BootCampMover();
        boot_mv.set_num_iterations(iters);
        TS_ASSERT( iters == boot_mv.get_num_iterations() )
    }

    void test_set_sfxn() {
        core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
        auto boot_mv = protocols::bootcamp::BootCampMover();
        boot_mv.set_sfxn(sfxn);
        TS_ASSERT( sfxn == boot_mv.get_sfxn() )
    }

    void test_xml_parse() {
        auto boot_mv = protocols::bootcamp::BootCampMover();
        //basic::resource_manager::ResourceManager rm;
        std::string xml_file = "<BootCampMover num_iterations=\"1025\"/>";
        std::string xml_file2 = "<BootCampMover scorefxn=\"testing123\"/>";
        //rm.read_resources_from_xml( "(unit test)", xml_file );
        utility::tag::TagCOP tag = tagptr_from_string(xml_file);
        utility::tag::TagCOP tag2 = tagptr_from_string(xml_file2);

        basic::datacache::DataMap data;
        core::scoring::ScoreFunctionOP sfxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
        //std::cout << "og name: " << sfxn->get_name() << std::endl;
        data.add( "scorefxns" , "testing123", sfxn );
        prime_Data( data );
        boot_mv.parse_my_tag(tag, data);
        //std::cout << boot_mv.get_num_iterations() << std::endl;
        TS_ASSERT( 1025 == boot_mv.get_num_iterations() )
        boot_mv.parse_my_tag(tag2, data);
        //std::cout << "og name: " << sfxn->get_name() << std::endl;
        //std::cout << "stored name: " << boot_mv.get_sfxn()->get_name() << std::endl;
        TS_ASSERT( sfxn == boot_mv.get_sfxn() )
        //TS_ASSERT_EQUALS( sfxn, boot_mv.get_sfxn() );
    }

    void test_bootcampmover_factory() {
        using namespace protocols::moves;

        std::string mover_name = "BootCampMover";
        MoverFactory * mover_factory = MoverFactory::get_instance();
        MoverOP my_mover = mover_factory->newMover( mover_name );

        protocols::bootcamp::BootCampMoverOP bcm_op =
                BootCampMoverOP( utility::pointer::dynamic_pointer_cast< protocols::bootcamp::BootCampMover > ( my_mover ) );

        //TS_ASSERT( identify_secondary_structure_spans( "" ).size() == 0 );
    }

};
