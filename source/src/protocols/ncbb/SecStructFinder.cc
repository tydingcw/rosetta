// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Project Headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/util.hh>

#include <core/kinematics/MoveMap.hh>

// Mover headers
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/ncbb/SecStructMinimizeMover.hh>
#include <protocols/ncbb/SecStructFinder.hh>
#include <protocols/ncbb/SecStructFinderCreator.hh>
#include <protocols/ncbb/util.hh>


//Basic headers
#include <basic/basic.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <string>
#include <sstream>

#include <core/conformation/Residue.hh> // AUTO IWYU For Pose::Residue Residue, Residue::ResidueType

// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace func;
using namespace constraints;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace core::pack::task;
using namespace core::id;


// tracer - used to replace cout
static basic::Tracer TR( "SecStructFinder" );

namespace protocols {
namespace ncbb {

SecStructFinder::SecStructFinder():
	residue_( "ALA" ),
	min_length_( 5 ),
	max_length_( 5 ),
	bin_size_( 10 ),
	dump_threshold_( -10000 ),
	dihedral_min_( -180 ),
	dihedral_max_(  180 ),
	dihedral_pattern_( expand_pattern_to_fit( "A", min_length_ ) ),
	alpha_beta_pattern_( expand_pattern_to_fit( "A", min_length_ ) ),
	min_everything_( false ),
	cart_( false ),
	constrain_( false ),
	dissimilarity_( 10 )
{
	Mover::type("SecStructFinder");

	score_fxn_ = get_score_function();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
}

SecStructFinder::SecStructFinder(
	std::string const & residue,
	core::Size min_length,
	core::Size max_length,
	Real bin_size,
	Real dissimilarity,
	Real dihedral_min,
	Real dihedral_max,
	Real dump_threshold,
	std::string dihedral_pattern,
	std::string alpha_beta_pattern,
	bool min_everything,
	bool cart,
	bool constrain
):
	residue_( residue ),
	min_length_( min_length ),
	//max_length_( max_length ),
	bin_size_( bin_size ),
	dump_threshold_( dump_threshold ),
	dihedral_min_( dihedral_min ),
	//dihedral_max_( ( dihedral_min > dihedral_max ) ? dihedral_min : dihedral_max ),
	//dihedral_pattern_( dihedral_pattern ),
	//alpha_beta_pattern_( expand_pattern_to_fit( alpha_beta_pattern, min_length_ ) ),
	min_everything_( min_everything ),
	cart_( cart ),
	constrain_( constrain ),
	dissimilarity_( ( bin_size > dissimilarity ) ? bin_size : dissimilarity )
{
	max_length_ = ( min_length_ > max_length ) ? min_length_ : max_length;
	dihedral_pattern_ = expand_pattern_to_fit( dihedral_pattern, min_length_ );
	alpha_beta_pattern_ = expand_pattern_to_fit( alpha_beta_pattern, min_length_ );
	dihedral_max_ = ( dihedral_min_ > dihedral_max ) ? dihedral_min_ : dihedral_max;

	Mover::type("SecStructFinder");

	score_fxn_ = get_score_function();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
}

std::string
SecStructFinder::alpha_to_beta( std::string alpha ) {
	if ( alpha == "ALA" ) {
		return "B3A";
	}
	if ( alpha == "CYS" ) {
		return "B3C";
	}
	if ( alpha == "ASP" ) {
		return "B3D";
	}
	if ( alpha == "GLU" ) {
		return "B3E";
	}
	if ( alpha == "PHE" ) {
		return "B3F";
	}
	if ( alpha == "GLY" ) {
		return "B3G";
	}
	if ( alpha == "HIS" ) {
		return "B3H";
	}
	if ( alpha == "ILE" ) {
		return "B3I";
	}
	if ( alpha == "LYS" ) {
		return "B3K";
	}
	if ( alpha == "LEU" ) {
		return "B3L";
	}
	if ( alpha == "MET" ) {
		return "B3M";
	}
	if ( alpha == "ASN" ) {
		return "B3N";
	}
	if ( alpha == "GLN" ) {
		return "B3Q";
	}
	if ( alpha == "ARG" ) {
		return "B3R";
	}
	if ( alpha == "SER" ) {
		return "B3S";
	}
	if ( alpha == "THR" ) {
		return "B3T";
	}
	if ( alpha == "VAL" ) {
		return "B3V";
	}
	if ( alpha == "TRP" ) {
		return "B3W";
	}
	if ( alpha == "TYR" ) {
		return "B3Y";
	}
	if ( alpha == "PRO" ) {
		return "B3P";
	}

	return "B3A";
}

std::string
SecStructFinder::expand_pattern_to_fit( std::string pattern, core::Size length ) {
	if ( length < pattern.length() ) {
		return pattern.substr( length );
	} else if ( length == pattern.length() ) {
		return pattern;
	} else {
		std::string expanded_pattern;
		while ( length > pattern.length() ) {
			expanded_pattern += pattern;
			length -= pattern.length();
		}
		expanded_pattern += pattern.substr( 0, length );
		return expanded_pattern;
	}
}

bool
SecStructFinder::uniq_refers_to_beta (
	char uniq
) {
	for ( core::Size i = 1; i <= dihedral_pattern_.length(); ++i ) {
		if ( dihedral_pattern_[i] == uniq ) {
			if ( alpha_beta_pattern_[i] == 'A' || alpha_beta_pattern_[i] == 'P' ) return false;
			else return true;
		}
	}
	return false;
}

void
SecStructFinder::initialize_pose( core::pose::Pose & pose ) {
	core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

	std::string name;
	std::string npatch;
	std::string cpatch;

	if ( alpha_beta_pattern_[ 0 ] == 'B' ) {
		name = alpha_to_beta( residue_ );
		npatch = ":AcetylatedNtermProteinFull";
	} else if ( alpha_beta_pattern_[ 0 ] == 'A' ) {
		name = residue_;
		npatch = ":AcetylatedNtermProteinFull";
	} else {
		name = "201";
		npatch = ":AcetylatedNtermPeptoidFull";
	}
	name += npatch;

	//restypes.push_back( residue_set_cap->name_map( name ) );
	pose.append_residue_by_jump( core::conformation::Residue( residue_set_cap->name_map( name ), true ), 1 );

	for ( core::Size resi = 1; resi < min_length_; ++resi ) {
		std::string name;
		std::string npatch;
		std::string cpatch;

		if ( alpha_beta_pattern_[ resi ] == 'B' ) {
			name = alpha_to_beta( residue_ );
			npatch = ":AcetylatedNtermProteinFull";
			cpatch = ":MethylatedCtermProteinFull";
		} else if ( alpha_beta_pattern_[ resi ] == 'A' ) {
			name = residue_;
			npatch = ":AcetylatedNtermProteinFull";
			cpatch = ":MethylatedCtermProteinFull";
		} else {
			name = "201";
			npatch = ":AcetylatedNtermPeptoidFull";
			cpatch = ":CtermPeptoidFull"; // can't constrain dihedral with peptoid final res
		}

		if ( resi == 0 ) {
			name += npatch;
		}

		// todo
		if ( resi == min_length_-1 ) {
			name += cpatch;
		}

		//restypes.push_back( residue_set_cap->name_map( name ) );
		pose.append_residue_by_bond( core::conformation::Residue( residue_set_cap->name_map( name ), true ), true );
	}
}

std::string
SecStructFinder::make_filename (
	core::Size number_dihedrals,
	utility::vector1< Real > dihedrals
) {
	std::stringstream filename;
	filename << residue_ << "_" << dihedral_pattern_ << "_" << min_length_;
	for ( core::Size i = 1; i <= number_dihedrals; ++i ) {
		filename << "_" << dihedrals[i];
	}
	filename << ".pdb";
	return filename.str();
}


bool
SecStructFinder::too_similar( core::Size i, core::Size j, utility::vector1< Real > dihedrals ) {
	bool similar = false;
	if ( alpha_beta_pattern_[i] == alpha_beta_pattern_[j] ) { //  don't compare dissimilars
		if ( alpha_beta_pattern_[i] == 'A' || alpha_beta_pattern_[i] == 'P' ) {
			if ( std::abs(basic::periodic_range( dihedrals[i] - dihedrals[j], 360 )) < dissimilarity_ && std::abs(basic::periodic_range( dihedrals[i+1] - dihedrals[j+1], 360 ) ) < dissimilarity_ ) {
				similar = true;
			}
		} else if ( alpha_beta_pattern_[i] == 'B' ) {
			if ( std::abs(basic::periodic_range( dihedrals[i] - dihedrals[j], 360 )) < dissimilarity_ && std::abs(basic::periodic_range( dihedrals[i+1] - dihedrals[j+1], 360 ))  < dissimilarity_ && std::abs(basic::periodic_range( dihedrals[i+2] - dihedrals[j+2], 360 ) )  < dissimilarity_ ) {
				similar = true;
			}
		}
	}
	return similar;
}

bool
increment_dihedrals(
	utility::vector1< Real > & dihedrals,
	Real const dihedral_min_,
	Real const dihedral_max_,
	Real const bin_size_
) {
	dihedrals[ 1 ] += bin_size_;
	core::Size p = 1;
	while ( dihedrals[ p ] >= dihedral_max_ ) {
		if ( p == dihedrals.size() ) return false;
		dihedrals[ p ] = dihedral_min_;
		dihedrals[ ++p ] += bin_size_;
	}
	return true;
}

void
SecStructFinder::show_current_dihedrals( core::Size number_dihedral_sets, utility::vector1< char > uniqs, utility::vector1< Real > dihedrals ) {
	TR << "Presently ";
	for ( core::Size i = 1; i <= number_dihedral_sets; ++i ) {
		TR << uniqs[i] << " ( ";
		core::Size take_from_here = give_dihedral_index( i, uniqs, dihedral_pattern_, alpha_beta_pattern_ );
		if ( uniq_refers_to_beta( uniqs[i] ) ) {
			TR << dihedrals[ take_from_here ] << ", " << dihedrals[ take_from_here+1 ] << ", " << dihedrals[ take_from_here+2 ] << " ), ";
		} else {
			TR << dihedrals[ take_from_here ] << ", " << dihedrals[ take_from_here+1 ] << " ), ";
		}
	}
}

void
SecStructFinder::apply( Pose & pose )
{
	core::Size number_dihedral_sets;
	utility::vector1< char > uniqs;
	count_uniq_char( dihedral_pattern_, number_dihedral_sets, uniqs );

	core::Size number_dihedrals = get_number_dihedrals( uniqs, dihedral_pattern_, alpha_beta_pattern_ );//2 * number_dihedral_sets; // plus one per beta AA!
	TR << "Investigating dihedral pattern " << dihedral_pattern_ << " with " << number_dihedral_sets << " uniques " << std::endl;

	initialize_pose( pose );

	kinematics::MoveMapOP min_mm( new kinematics::MoveMap );
	min_mm->set_bb( true );
	min_mm->set_chi( true );

	protocols::minimization_packing::MinMover minmover( min_mm, score_fxn_, "lbfgs_armijo_nonmonotone", 0.0001, true );
	minmover.cartesian( cart_ );
	utility::vector1< Pose > poses_to_min;
	utility::vector1< utility::vector1< Real > > dihedrals_to_min;
	utility::vector1< std::string > minned_filenames;

	// now we have our pose
	utility::vector1< Real > dihedrals( number_dihedrals+1, dihedral_min_ );

	// while loop idiom for looping through every layer of a vector of indices
	do {

		// skip if the dihedrals for A are the same as B, etc.
		// break after first example found
		bool skip = false;
		for ( core::Size i = 1; i <= number_dihedral_sets; i += 2 ) {
			for ( core::Size j = i+1; j <= number_dihedral_sets; j++ ) {
				skip = too_similar( i, j, dihedrals );
				if ( skip ) break;
			}
			if ( skip ) break;
		}

		if ( !skip ) {

			TR << "Trying dihedral DOF vector: ";
			for ( core::Size ii = 1; ii <= number_dihedrals; ++ii ) {
				TR << dihedrals[ii] << ", ";
			} TR << std::endl;

			for ( core::Size resi = 1; resi <= pose.size(); ++resi ) {
				core::Size take_from_here = 1;
				for ( core::Size i = 1; i <= number_dihedral_sets; ++i ) {
					if ( dihedral_pattern_[ resi-1 ] == uniqs[ i ] ) { // -1 to sync indexing
						take_from_here = give_dihedral_index( i, uniqs, dihedral_pattern_, alpha_beta_pattern_  );
					}
				}

				//TR << "For res " << resi << " ( " << pose.residue_type( resi ).name() << " )... " << std::endl;
				if ( pose.residue( resi ).type().is_beta_aa() ) {
					//TR << "Assigning DOF vector items " << take_from_here << " to " << take_from_here+2 << " to this BAA" << std::endl;
					pose.set_torsion( TorsionID( resi, BB, 1 ), dihedrals[ take_from_here ] );
					pose.set_torsion( TorsionID( resi, BB, 2 ), dihedrals[ take_from_here+1 ] );
					pose.set_torsion( TorsionID( resi, BB, 3 ), dihedrals[ take_from_here+2 ] );
					pose.set_torsion( TorsionID( resi, BB, 4 ), 180 );
				} else {
					//TR << "Assigning DOF vector items " << take_from_here << " to " << take_from_here+1 << " to this AAA" << std::endl;
					pose.set_phi( resi, dihedrals[ take_from_here ] ); //phi_vec[ vec_index ] );
					pose.set_psi( resi, dihedrals[ take_from_here+1 ] ); //psi_vec[ vec_index ] );
					pose.set_omega( resi, 180 );
				}
			}

			show_current_dihedrals( number_dihedral_sets, uniqs, dihedrals );

			Real score = ( *score_fxn_ ) ( pose );
			TR << "score is " << score << std::endl;
			score_fxn_->show(pose);
			if ( min_everything_ ) {

				if ( score_fxn_->has_zero_weight( core::scoring::dihedral_constraint ) ) {
					score_fxn_->set_weight( core::scoring::dihedral_constraint, 1.0 );
				}

				Pose minpose = pose;

				SecStructMinimizeMoverOP ssmm( new SecStructMinimizeMover() );
				ssmm->set_constrain( constrain_ );
				ssmm->set_scorefunction( score_fxn_ );
				ssmm->set_alpha_beta_pattern( alpha_beta_pattern_ );
				ssmm->set_dihedral_pattern( dihedral_pattern_ );
				ssmm->set_dihedrals( dihedrals );
				ssmm->apply( minpose );

				score_fxn_->set_weight( core::scoring::dihedral_constraint, 0.0 );
				Real score = ( ( *score_fxn_ ) ( minpose ) );
				score_fxn_->show(minpose);
				TR << " and minned is " << score << std::endl;
				if ( score <= dump_threshold_ ) {
					TR << "Thus, dumping ";
					std::string filename = make_filename( number_dihedrals, dihedrals );

					minpose.dump_scored_pdb( filename.c_str(), ( *score_fxn_ ) );
					pose.dump_scored_pdb( "fuck.pdb", ( *score_fxn_ ) );

					if ( minpose.residue_type( 1 ).is_beta_aa() ) {
						TR << " ( " << minpose.torsion( TorsionID( 1, BB, 1 ) ) << ", "<< minpose.torsion( TorsionID( 1, BB, 2 ) ) << ", " << minpose.torsion( TorsionID( 1, BB, 3 ) ) << " ), ";
					} else {
						TR << " ( " << minpose.torsion( TorsionID( 1, BB, 1 ) ) << ", " << minpose.psi( 1 ) << " ), ";
					}

					for ( core::Size resi = 2; resi < pose.size(); ++resi ) {
						if ( minpose.residue_type( resi ).is_beta_aa() ) {
							TR << " ( " << minpose.torsion( TorsionID( resi, BB, 1 ) ) << ", "<< minpose.torsion( TorsionID( resi, BB, 2 ) ) << ", " << minpose.torsion( TorsionID( resi, BB, 3 ) ) << " ), ";
						} else {
							TR << " ( " << minpose.phi( resi ) << ", " << minpose.psi( resi ) << " ), ";
						}
					}

					if ( minpose.residue_type( pose.size() ).is_beta_aa() ) {
						TR << " ( " << minpose.torsion( TorsionID(pose.size(), BB, 1 ) ) << ", "<< minpose.torsion( TorsionID( pose.size(), BB, 2 ) ) << ", "
							<< minpose.torsion( TorsionID( pose.size(), BB, 3 ) ) << " ), ";
					} else {
						TR << " ( " << minpose.phi( pose.size() ) << ", " << minpose.torsion( TorsionID( pose.size(), BB, 2 ) ) << " ), " << std::endl;
					}
				}
			} else {
				if ( score <= dump_threshold_ ) {
					TR << " thus, dumping ";
					std::string filename = make_filename( number_dihedrals, dihedrals );

					pose.dump_scored_pdb( filename.c_str(), ( *score_fxn_ ) );
					TR << filename.c_str() << " " << score << std::endl;
					poses_to_min.push_back( pose );
					utility::vector1< Real > dihedrals_for_minimization;
					for ( core::Size index = 1; index <= number_dihedrals; ++index ) {
						dihedrals_for_minimization.push_back( dihedrals[ index ] );
					}
					dihedrals_to_min.push_back( dihedrals_for_minimization );
					minned_filenames.push_back( filename );
				}
			}
		}
	} while ( increment_dihedrals( dihedrals, dihedral_min_, dihedral_max_, bin_size_ ) );


	// now min all our poses to min
	if ( poses_to_min.size() == 0 ) return;

	for ( core::Size ii = 1, sz = poses_to_min.size(); ii < sz; ++ii ) {

		Pose minpose = poses_to_min[ii];

		SecStructMinimizeMoverOP ssmm( new SecStructMinimizeMover() );
		ssmm->set_constrain( constrain_ );
		ssmm->set_scorefunction( score_fxn_ );
		ssmm->set_alpha_beta_pattern( alpha_beta_pattern_ );
		ssmm->set_dihedral_pattern( dihedral_pattern_ );
		ssmm->set_dihedrals( dihedrals );
		ssmm->apply( minpose );

		score_fxn_->set_weight( core::scoring::dihedral_constraint, 0.0 );
		Real score = ( ( *score_fxn_ ) ( minpose ) );
		TR << "for stored pose with dihedrals ";
		for ( core::Size resi = 1, n = poses_to_min[ii].size(); resi <= n; ++resi ) {
			if ( alpha_beta_pattern_[ resi-1 ] == 'A' || alpha_beta_pattern_[ resi-1 ] == 'P' ) {
				TR << "( " << poses_to_min[ii].phi( resi ) << ", " << poses_to_min[ii].psi( resi ) << " ), ";
			} else {
				Real phi = poses_to_min[ii].torsion( TorsionID( resi, BB, 1 )  );
				Real tht = poses_to_min[ii].torsion( TorsionID( resi, BB, 2 )  );
				Real psi = poses_to_min[ii].torsion( TorsionID( resi, BB, 3 )  );
				TR << "( " << phi << ", " << tht << ", " << psi << " ), ";
			}
		}
		TR << std::endl << " now ";
		for ( core::Size resi = 1, n = minpose.size(); resi <= n; ++resi ) {
			if ( alpha_beta_pattern_[ resi-1 ] == 'A' || alpha_beta_pattern_[ resi-1 ] == 'P' ) {
				TR << "( " << minpose.phi( resi ) << ", " << minpose.psi( resi ) << " ), ";
			} else {
				Real phi = minpose.torsion( TorsionID( resi, BB, 1 )  );
				Real tht = minpose.torsion( TorsionID( resi, BB, 2 )  );
				Real psi = minpose.torsion( TorsionID( resi, BB, 3 )  );
				TR << "( " << phi << ", " << tht << ", " << psi << " ), ";
			}
		}

		TR << "minned score is " << score << std::endl;
		minned_filenames[ii] = "minned_" + minned_filenames[ii];
		minpose.dump_scored_pdb( minned_filenames[ii].c_str(), ( * score_fxn_ ) );
	}
}//apply

protocols::moves::MoverOP SecStructFinder::clone() const {
	return utility::pointer::make_shared< SecStructFinder >(
		residue_,
		min_length_,
		max_length_,
		bin_size_,
		dissimilarity_,
		dihedral_min_,
		dihedral_max_,
		dump_threshold_,
		dihedral_pattern_,
		alpha_beta_pattern_,
		min_everything_,
		cart_ );
}

void SecStructFinder::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
) {

	residue_ = tag->getOption<std::string>("residue", "ALA" );
	min_length_ = tag->getOption<core::Size>("min_length", 5 );
	max_length_ = tag->getOption<core::Size>("max_length", 5 );

	if ( max_length_ < min_length_ ) {
		TR.Warning << "Provided max_length was less than min_length, setting max equal to min" << std::endl;
		max_length_ = min_length_;
	}

	bin_size_ = tag->getOption<Real>( "bin_size", 10.0 );
	dissimilarity_ = tag->getOption<Real>( "dissimilarity", 10.0 );

	if ( dissimilarity_ < bin_size_ ) {
		TR.Warning << "Provided a dissimilarity score less than the bin size, which is meaningless. Setting equal to bin size." << std::endl;
		dissimilarity_ = bin_size_;
	}

	dump_threshold_ = tag->getOption<Real>( "dump_threshold", 10.0 );
	dihedral_pattern_ = tag->getOption<std::string>( "dihedral_pattern", "A" );

	TR << "Expanding dihedral pattern to fit the length provided: ";
	dihedral_pattern_ = expand_pattern_to_fit( dihedral_pattern_, min_length_ );
	TR << dihedral_pattern_ << std::endl;

	alpha_beta_pattern_ = tag->getOption<std::string>("alpha_beta_pattern", "A" );

	TR << "Expanding alpha/beta pattern to fit the length provided: ";
	alpha_beta_pattern_ = expand_pattern_to_fit( alpha_beta_pattern_, min_length_ );
	TR << alpha_beta_pattern_ << std::endl;

	min_everything_ = tag->getOption<bool>( "min_everything", false );
	cart_ = tag->getOption<bool>("cart", false );
}

// MoverCreator
moves::MoverOP SecStructFinderCreator::create_mover() const {
	return utility::pointer::make_shared< SecStructFinder >();
}

std::string SecStructFinderCreator::keyname() const {
	return SecStructFinderCreator::mover_name();
}

std::string SecStructFinderCreator::mover_name(){
	return "SecStructFinder";
}

}
}
