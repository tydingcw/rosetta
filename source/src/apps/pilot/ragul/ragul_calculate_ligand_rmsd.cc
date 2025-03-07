// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Ragul Gowthaman

//GPU enabling is not default
//To test how many threads are fastest for your computer,
//use -gpu:threads 1024 (or other number) on the command line

#include <iostream>
#include <string>
#include <cmath>

// Protocol Headers
#include <devel/init.hh>
#include <basic/options/option_macros.hh>

// Utility Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>



#include <core/import_pose/import_pose.hh>
#include <utility/options/StringOption.hh>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace basic::options::OptionKeys;

OPT_KEY( String, input_ligand )
OPT_KEY( String, reference_ligand )

int main( int argc, char * argv [] ) {

	try{

		NEW_OPT( input_ligand, "ligand file name", "input_ligand.pdb" );
		NEW_OPT( reference_ligand, "reference ligand file name", "reference_ligand.pdb" );

		devel::init(argc, argv);

		//setup reference ligand
		pose::Pose ref_pose;
		std::string const ref_ligand = option[ reference_ligand ];
		core::import_pose::pose_from_file( ref_pose, ref_ligand , core::import_pose::PDB_file);
		core::Size ref_res_num = 0;
		for ( int j = 1, resnum = ref_pose.size(); j <= resnum; ++j ) {
			if ( !ref_pose.residue(j).is_protein() ) {
				ref_res_num = j;
				break;
			}
		}
		if ( ref_res_num == 0 ) {
			std::cout<<"Error, no reference ligand for RMSD calculation" << std::endl;
			exit(1);
		}
		core::conformation::Residue const & pose1_rsd = ref_pose.conformation().residue(ref_res_num);

		//setup input ligand
		pose::Pose inp_pose;
		std::string const inp_ligand = option[ input_ligand ];
		core::import_pose::pose_from_file( inp_pose, inp_ligand , core::import_pose::PDB_file);
		core::Size inp_res_num = 0;
		for ( int j = 1, resnum = inp_pose.size(); j <= resnum; ++j ) {
			if ( !inp_pose.residue(j).is_protein() ) {
				inp_res_num = j;
				break;
			}
		}
		if ( inp_res_num == 0 ) {
			std::cout<<"Error, no input ligand for RMSD calculation" << std::endl;
			exit(1);
		}
		core::conformation::Residue const & pose2_rsd = inp_pose.conformation().residue(inp_res_num);


		//CALCULATE RMSD
		//IMPORTANT:this rmsd calculation does not consider symmetry
		//assert input & reference ligand should have same num of atoms

		if ( !( pose2_rsd.nheavyatoms() == pose1_rsd.nheavyatoms() ) ) {
			std::cout<<"Error, input & reference ligand should have same no. of atoms in the same order" << std::endl;
			exit(1);
		} else {
			core::Real rmsd, dist_sum(0.);
			Size j = 1;
			for ( Size i = 1; i <= pose1_rsd.nheavyatoms(); ++i ) {
				if ( j <= pose2_rsd.nheavyatoms() ) {
					core::Real x_dist =  ( (pose1_rsd.atom(i).xyz()(1) - pose2_rsd.atom(j).xyz()(1)) * (pose1_rsd.atom(i).xyz()(1) - pose2_rsd.atom(j).xyz()(1)) );
					core::Real y_dist =  ( (pose1_rsd.atom(i).xyz()(2) - pose2_rsd.atom(j).xyz()(2)) * (pose1_rsd.atom(i).xyz()(2) - pose2_rsd.atom(j).xyz()(2)) );
					core::Real z_dist =  ( (pose1_rsd.atom(i).xyz()(3) - pose2_rsd.atom(j).xyz()(3)) * (pose1_rsd.atom(i).xyz()(3) - pose2_rsd.atom(j).xyz()(3)) );
					dist_sum += x_dist + y_dist + z_dist;
				}
				++j;
			}
			rmsd = sqrt(dist_sum/pose1_rsd.nheavyatoms());
			std::cout<<"RMSD:"<<"\t"<<inp_ligand<<"\t"<<ref_ligand<<"\t"<<rmsd<<std::endl;
		}
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
