// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk

#include <iostream>

// Protocol Headers

// Core Headers
#include <core/conformation/Conformation.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/options/util.hh>
//#include <protocols/pockets/PocketConstraint.hh>
#include <protocols/pockets/PocketGrid.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

// Utility Headers

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>

#include <fstream> // AUTO IWYU For filebuf


using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace basic::options::OptionKeys;

OPT_KEY( Integer, num_angles )

/// General testing code
int
main( int argc, char * argv [] )
{
	try {

		NEW_OPT ( num_angles, "Number of different pose angles to measure score at", 1);

		char chain = '\0';
		devel::init(argc, argv);
		int angles = option[ num_angles ];
		if ( angles <1 ) {
			fprintf (stderr, "Error: invalid number of angles.  Must be greather than 0\n");
			return -1;
		}

		std::string const output_tag = option[ OptionKeys::out::output_tag ]();
		pose::Pose input_pose;

		//read in pdb file from command line
		std::string const input_pdb_name ( basic::options::start_file() );
		core::import_pose::pose_from_file( input_pose, input_pdb_name , core::import_pose::PDB_file);

		//scoring::ScoreFunctionOP scorefxn( get_score_function() );

		// report starting score
		//(*scorefxn)(input_pose);
		//core::Real const starting_total_score = input_pose.energies().total_energies()[ total_score ];
		//std::cout << "Total score at start without constraint is: " << starting_total_score << std::endl;
		//core::Real const starting_pocket_score = input_pose.energies().total_energies()[ pocket_constraint ];
		//std::cout << "Constraint score (unweighted) at start without constraint is: " << starting_pocket_score << std::endl;

		// create pocket constraint
		//core::scoring::constraints::PocketConstraintOP pcons ( new core::scoring::constraints::PocketConstraint(input_pose) );
		//scorefxn->set_weight( core::scoring::pocket_constraint, 1 );

		std::filebuf fb,fb2;
		std::stringstream filename,filename2;
		if ( !option[ OptionKeys::out::output_tag ]().empty() ) {
			filename<<option[ OptionKeys::out::output_tag ]()<<".pscore";
			filename2<<option[ OptionKeys::out::output_tag ]()<<".lpscore";
		} else {
			filename<<basic::options::start_file()<<".pscore";
			filename2<<basic::options::start_file()<<".lpscore";
		}
		fb.open (filename.str().c_str(),std::ios::out);
		fb2.open (filename2.str().c_str(),std::ios::out);
		std::ostream os(&fb);
		std::ostream os2(&fb2);

		//Loop over all residues and get score
		for ( int j = 1, resnum = input_pose.size(); j <= resnum; ++j ) {
			if ( chain == '\0' ) {
				chain = input_pose.pdb_info()->chain(j);
			}
			if ( chain != input_pose.pdb_info()->chain(j) ) break;
			//std::cout << input_pose.pdb_info()->chain(j)<<input_pose.pdb_info()->number(j)<< std::endl;
			//pcons->set_target_res(input_pose, j);
			//input_pose.add_constraint( pcons );
			//(*scorefxn)(input_pose);
			//core::Real constraint_pocket_score = input_pose.energies().total_energies()[ pocket_constraint ];
			core::Real constraint_pocket_score = 0;
			core::Real largest_pocket_score = 0;
			for ( int i=0; i<angles; ++i ) {
				if ( i>0 ) {
					core::Real x,y,z;
					x = (int) (numeric::random::uniform() *89 +1);
					y = (int) (numeric::random::uniform() *89 +1);
					z = (int) (numeric::random::uniform() *89 +1);
					numeric::xyzMatrix<core::Real> x_rot_mat( numeric::x_rotation_matrix_degrees(x) );
					numeric::xyzMatrix<core::Real> y_rot_mat( numeric::y_rotation_matrix_degrees(y) );
					numeric::xyzMatrix<core::Real> z_rot_mat( numeric::z_rotation_matrix_degrees(z) );
					core::Vector v(0,0,0);
					input_pose.apply_transform_Rx_plus_v(x_rot_mat, v);
					input_pose.apply_transform_Rx_plus_v(y_rot_mat, v);
					input_pose.apply_transform_Rx_plus_v(z_rot_mat, v);
				}

				protocols::pockets::PocketGrid pg( input_pose.conformation().residue(j) );
				if ( pg.autoexpanding_pocket_eval( input_pose.conformation().residue(j), input_pose ) ) {
					constraint_pocket_score += pg.netTargetPocketVolume();
					largest_pocket_score += pg.largestTargetPocketVolume();
				}
			}

			constraint_pocket_score /= angles;
			largest_pocket_score /= angles;
			os << input_pose.pdb_info()->chain(j)<<input_pose.pdb_info()->number(j)<<"\t" << constraint_pocket_score << std::endl;
			os2 << input_pose.pdb_info()->chain(j)<<input_pose.pdb_info()->number(j)<<"\t" << largest_pocket_score << std::endl;
			std::cout << input_pose.pdb_info()->chain(j)<<input_pose.pdb_info()->number(j)<<"\t" << constraint_pocket_score << std::endl;
			std::cout << input_pose.pdb_info()->chain(j)<<input_pose.pdb_info()->number(j)<<"\tLargest: " << largest_pocket_score << std::endl;
			//input_pose.remove_constraints( );
		}

		fb.close();
		fb2.close();
		/*
		//input_pose.remove_constraint( pcons );
		pcons->set_target_res_pdb(input_pose, "53");
		input_pose.add_constraint( pcons );
		(*scorefxn)(input_pose);
		constraint_pocket_score = input_pose.energies().total_energies()[ pocket_constraint ];
		std::cout << "Constraint score (unweighted) at start with constraint is: " << constraint_pocket_score << std::endl;
		input_pose.remove_constraints( );
		//input_pose.remove_constraint( pcons );
		pcons->set_target_res_pdb(input_pose, "54");
		input_pose.add_constraint( pcons );
		(*scorefxn)(input_pose);
		constraint_pocket_score = input_pose.energies().total_energies()[ pocket_constraint ];
		std::cout << "Constraint score (unweighted) at start with constraint is: " << constraint_pocket_score << std::endl;
		input_pose.remove_constraints( );
		//input_pose.remove_constraint( pcons );
		pcons->set_target_res_pdb(input_pose, "55");
		input_pose.add_constraint( pcons );
		(*scorefxn)(input_pose);
		constraint_pocket_score = input_pose.energies().total_energies()[ pocket_constraint ];
		std::cout << "Constraint score (unweighted) at start with constraint is: " << constraint_pocket_score << std::endl;

		//core::scoring::constraints::PocketConstraint pcons( input_pose );
		//input_pose.add_constraint( pcons );

		// rescore, report new score
		//(*scorefxn)(input_pose);
		//core::Real const constraint_total_score = input_pose.energies().total_energies()[ total_score ];
		//std::cout << "Total score at start with constraint is: " << constraint_total_score << std::endl;
		//core::Real const constraint_pocket_score = input_pose.energies().total_energies()[ pocket_constraint ];
		//std::cout << "Constraint score (unweighted) at start with constraint is: " << constraint_pocket_score << std::endl;

		*/

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}


