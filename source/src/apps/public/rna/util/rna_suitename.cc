// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <basic/options/option.hh>
#include <protocols/viewer/viewers.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/rna/RNA_SuiteName.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <devel/init.hh>
#include <ObjexxFCL/format.hh>
#include <protocols/rna/denovo/util.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL::format;

///////////////////////////////////////////////////////////////////////////////
void
rna_suitename()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::import_pose::pose_stream;
	using namespace core::chemical;
	using namespace core::pose::rna;
	using namespace protocols::rna::denovo;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// input stream
	PoseInputStreamOP input;
	input = utility::pointer::make_shared< PDBPoseInputStream >( option[ in::file::s ]() );


	pose::Pose pose;
	Size i( 0 );
	RNA_SuiteName suitename;

	while ( input->has_another_pose() ) {
		input->fill_pose( pose, *rsd_set );
		i++;

		//  protocols::rna::denovo::ensure_phosphate_nomenclature_matches_mini( pose );
		core::pose::rna::figure_out_reasonable_rna_fold_tree( pose );
		core::pose::rna::virtualize_5prime_phosphates( pose ); // should we have this on by deafult?

		std::cout << "-----Pose " << i << ": " << tag_from_pose( pose ) << "-----" << std::endl;
		for ( Size j = 1; j <= pose.size(); ++j ) {
			if ( !pose.residue( j ).is_RNA() ) continue;
			RNA_SuiteAssignment assignment = suitename.assign(pose, j + 1 ); // note funny j+1 convention in suitename.assign
			std::cout << "Residue " << I( 3, j ) << ' ' <<
				pose.pdb_info()->chain(j) << ':' << I(4, pose.pdb_info()->number(j)) <<
				' ' << pose.residue(j).name3() << ' ' <<
				assignment.name << ' ' << std::setprecision(3) <<assignment.suiteness << std::endl;
		}
	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	rna_suitename();

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -s <pdb file> " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		option.add_relevant( in::file::s );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
}
