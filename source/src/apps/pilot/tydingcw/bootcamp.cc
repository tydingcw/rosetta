// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <iostream>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/random/random.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include  <protocols/bootcamp/fold_tree_from_ss.hh>

int main( int argc, char ** argv ) {

	using namespace protocols::moves;

	std::cout << "Hello World!" << std::endl;

	devel::init( argc, argv );
	utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
	if ( filenames.size() > 0 ) {
		std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
	} else {
		std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
		return 1;
	}

	core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	core::Real score = sfxn->score( *mypose );
	std::cout << "Score: " << score << std::endl;

    //modify foldtree
    auto my_tree = protocols::bootcamp::fold_tree_from_ss( *mypose );
    mypose->fold_tree(my_tree);

    //chainbreak terms
    sfxn->set_weight(core::scoring::linear_chainbreak, 1);
    core::pose::correctly_add_cutpoint_variants(*mypose);

	core::Size n_resi = mypose->size();

	//core::Real gauss = numeric::random::gaussian();
	//std::cout << gauss << std::endl;

	core::Real uni = numeric::random::uniform();
	//std::cout << uni << std::endl;

	//core::Size randres = n_resi * uni + 1;
	//std::cout << randres << " randres"  << std::endl;


	MonteCarlo monte_carlo = MonteCarlo(*mypose, *sfxn, 0.6);

	protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( *mypose, true, 0 );
	the_observer->pymol().apply( *mypose);

	core::kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_chi( true );

	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );

	core::optimization::AtomTreeMinimizer atm;

	core::pose::Pose copy_pose;

    bool accept;
    core::Real accept_count = 0;
    core::Real accept_ratio;
    core::Real sum_scores = 0.0;
    core::Real avg_score;

	for (int i = 1; i <= 100; i++){

	    core::Size randres = n_resi * uni + 1;
	    core::Real pert1 = numeric::random::gaussian();
	    core::Real pert2 = numeric::random::gaussian();

	    core::Real orig_phi = mypose->phi( randres );
	    core::Real orig_psi = mypose->psi( randres );
	    mypose->set_phi( randres, orig_phi + pert1 );
	    mypose->set_psi( randres, orig_psi + pert2 );

	    core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
	    repack_task->restrict_to_repacking();
	    core::pack::pack_rotamers( *mypose, *sfxn, repack_task );

	    //atm.run( *mypose, mm, *sfxn, min_opts );

	    copy_pose = *mypose;
	    atm.run( copy_pose, mm, *sfxn, min_opts );
	    *mypose = copy_pose;

        //Record Acceptance ratio and score
	    accept = monte_carlo.boltzmann(*mypose);
        if (accept) {
            accept_count++;
        }
        score = sfxn->score( *mypose );
        sum_scores += score;
        //Acceptance Ratio at Step 100: 0.81
        //Average Score at Step 100: -246.001
        if (i%100 == 0) {
            accept_ratio = accept_count/i;
            std::cout << "Acceptance Ratio at Step " << i << ": " << accept_ratio << std::endl;
            avg_score = sum_scores/i;
            std::cout << "Average Score at Step " << i << ": " << avg_score << std::endl;
        }
		
	}

	return 0;
}

