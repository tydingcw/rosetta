// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/calc_taskop_movers/DesignRepackMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Eva-Maria Strauch (evas01@u.washington.edu)
#include <protocols/calc_taskop_movers/DesignRepackMover.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <protocols/scoring/Interface.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/PackerTask.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/make_symmetric_task.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/FalseResidueSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/util.hh>

#include <protocols/moves/Mover.hh>
//#include <protocols/moves/ResidueMover.hh>


#include <basic/options/option.hh>

#include <utility/tag/Tag.hh>

// Utility Headers

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>

// Unit Headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>
#include <core/pose/selection.hh>


// C++ headers
#include <string>

#include <core/conformation/Residue.hh> // AUTO IWYU For Pose::Residue Residue, Residue::ResidueType
#include <utility/tag/XMLSchemaGeneration.hh> // AUTO IWYU For operator+, XMLSchemaAttribute, xsct_...
#include <core/pack/task/ResidueLevelTask.hh> // AUTO IWYU For ResidueLevelTask

using namespace core;

// option key includes

#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <protocols/task_operations/ProteinInterfaceDesignOperation.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>


using namespace core::scoring;

static basic::Tracer TR( "protocols.simple_moves.DesignRepackMover" );


namespace protocols {
namespace calc_taskop_movers {

using namespace protocols::moves;

DesignRepackMover::DesignRepackMover() : protocols::moves::Mover( "DesignRepackMover" ),
	repack_partner1_( false ), repack_partner2_( false ), design_partner1_( false ), design_partner2_( false ), min_rb_set_( false ), interface_distance_cutoff_(8.0 ), repack_non_ala_( true ), optimize_foldtree_( true ), automatic_repacking_definition_( true ), use_preset_task_( false ), symmetry_( false )
{
	allowed_aas_.resize( core::chemical::num_canonical_aas, true );
	prevent_repacking_.clear();
	restrict_to_repacking_.clear();
}

DesignRepackMover::DesignRepackMover( std::string const & name ) : protocols::moves::Mover( name ), repack_partner1_( false ), repack_partner2_( false ), design_partner1_( false ), design_partner2_( false ), min_rb_set_( false ), interface_distance_cutoff_( 8.0 ), repack_non_ala_( true ), optimize_foldtree_( true ), automatic_repacking_definition_( true ), use_preset_task_( false ), symmetry_( false )
{
	allowed_aas_.resize( core::chemical::num_canonical_aas, true );
}

std::string
DesignRepackMover::get_name() const {
	return "DesignRepackMover";
}

void
DesignRepackMover::clear_task(){
	task_ = nullptr;
}

void
DesignRepackMover::min_sc( core::select::residue_selector::ResidueSelectorOP min_sc ) { min_sc_ = min_sc; }

utility::vector1< bool >
DesignRepackMover::min_sc( core::pose::Pose const & pose ) const {
	runtime_assert( min_sc_ != nullptr );
	return min_sc_->apply( pose );
}

bool
DesignRepackMover::min_sc_set() const {
	return min_sc_ != nullptr;
}

void
DesignRepackMover::min_bb( core::select::residue_selector::ResidueSelectorOP min_bb ) { min_bb_ = min_bb; }

utility::vector1< bool >
DesignRepackMover::min_bb( core::pose::Pose const & pose ) const {
	runtime_assert( min_bb_ != nullptr );
	return min_bb_->apply( pose );
}

bool
DesignRepackMover::min_bb_set() const {
	return min_bb_ != nullptr;
}

void
DesignRepackMover::task_factory( core::pack::task::TaskFactoryOP p ) {
	task_factory_ = p;
}

void
DesignRepackMover::setup_packer_and_movemap( core::pose::Pose const & in_pose )
{
	core::Size const rb_jump( 1 );

	core::pose::Pose pose( in_pose ); // to maintain constness of input pose
	curr_min_bb_.clear(); curr_min_sc_.clear(); curr_min_rb_.clear();
	curr_min_bb_.resize( pose.size(), false );
	curr_min_sc_.resize( pose.size(), false );
	curr_min_rb_.resize( pose.num_jump(), true );
	if ( min_sc_set() ) {
		curr_min_sc_ = min_sc( pose );
	}
	if ( min_bb_set() ) {
		curr_min_bb_ = min_bb( pose );
	}
	if ( min_rb_set() ) curr_min_rb_ = min_rb_;

	using ObjexxFCL::FArray1D_bool;
	FArray1D_bool partner1( pose.size() );
	if ( !symmetry_ ) {
		pose.fold_tree().partition_by_jump( rb_jump, partner1 ); // partner1 is true for all residues in partner1; false o/w
	}

	protocols::scoring::Interface interface_obj(rb_jump);
	pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`
	interface_obj.distance( interface_distance_cutoff_ );
	interface_obj.calculate( pose );

	if ( !task_ || !use_preset_task_ ) {
		/// This allows presetting the task
		/// note that the operations on desiganble positions are of logical intersection type and only restrict the identities
		/// allowed at the positions, so if more rotamers and less identities are set before calling this setup, they
		/// will be respected
		if ( task_factory_ ) {
			task_ = task_factory_->create_task_and_apply_taskoperations( pose );
		} else {
			task_ = core::pack::task::TaskFactory::create_packer_task( pose );
		}
		task_->initialize_from_command_line();
	}

	protocols::task_operations::ProteinInterfaceDesignOperation pido;
	pido.repack_chain1( repack_partner1_ );
	pido.repack_chain2( repack_partner2_ );
	pido.design_chain1( design_partner1_ );
	pido.design_chain2( design_partner2_ );
	pido.interface_distance_cutoff( interface_distance_cutoff_ );
	pido.apply( pose, *task_ );

	// Restrict to allowable amino acids
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		core::pack::task::operation::RestrictAbsentCanonicalAAS racaas( i, allowed_aas_ );
		racaas.apply( pose, *task_ );
	}

	utility::vector1< core::Size > target_res = target_residues( pose );

	// Packertask defaults to designing all residues, so this paragraph only has to restrict
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		if ( interface_obj.is_interface( i ) ) { // in interface
			if ( (partner1( i ) && repack_partner1_ ) || (!partner1( i ) && repack_partner2_) ) { //redesign/repack interface
				core::Size const restype( pose.residue(i).aa() );
				// Check for pro, gly, and disulfide bonded cysteines
				if ( ( !repack_non_ala_ && restype != chemical::aa_ala ) || restype == chemical::aa_pro || restype == chemical::aa_gly || pose.residue(i).type().is_disulfide_bonded() ) {
					//     if( automatic_repacking_definition_ ) task_->nonconst_residue_task(i).prevent_repacking();
					if ( ! ( pose.residue( i ).type().is_disulfide_bonded() ) ) {
						if ( !min_sc_set() ) curr_min_sc_[ i ] = true;
					}
					if ( !min_bb_set() ) curr_min_bb_[ i ] = true;
				} else { // non-PRO/GLY/CYD
					if ( target_res.size() == 0 ) { // target residues undefined => repack the entire interface
						//      if( design_ && automatic_repacking_definition_ )
						//       task_->nonconst_residue_task(i).restrict_absent_canonical_aas( allowed_aas_ );
						//      else if( !design_ && automatic_repacking_definition_ )
						//       task_->nonconst_residue_task( i ).restrict_to_repacking();
						if ( !min_sc_set() ) curr_min_sc_[ i ] = true;
						if ( !min_bb_set() ) curr_min_bb_[ i ] = true;
					} else { //target residues defined
						core::conformation::Residue const resi( pose.residue( i ) );

						for ( utility::vector1< core::Size >::const_iterator target_it = target_res.begin();
								target_it!=target_res.end(); ++target_it ) {
							core::conformation::Residue const res_target( pose.residue( *target_it ) );

							Real const distance( resi.xyz( resi.nbr_atom() ).distance( res_target.xyz( res_target.nbr_atom() ) ) );
							if ( distance <= interface_distance_cutoff_ ) {
								//        if( design_ && automatic_repacking_definition_ )
								//         task_->nonconst_residue_task( i ).restrict_absent_canonical_aas( allowed_aas_ );
								//        else if( !design_ && automatic_repacking_definition_ )
								//         task_->nonconst_residue_task( i ).restrict_to_repacking();
								if ( !min_sc_set() ) curr_min_sc_[ i ] = true;
								if ( !min_bb_set() ) curr_min_bb_[ i ] = true;

								break;
							} else { //within distance cutoff
								//         if( automatic_repacking_definition_ )
								//         task_->nonconst_residue_task(i).prevent_repacking();
								if ( !min_sc_set() ) curr_min_sc_[ i ] = false;
								if ( !min_bb_set() ) curr_min_bb_[ i ] = true;
							}
						} // target_res iterator
						// reach this point only if target residues are defined, and current residue is at the interface, but not close to target
						// and for some reason was not taken care of before the break statement above
						//      if( automatic_repacking_definition_ )
						//       task_->nonconst_residue_task( i ).prevent_repacking();
						if ( !min_sc_set() ) curr_min_sc_[ i ] = false;
						if ( !min_bb_set() ) curr_min_bb_[ i ] = true;
					} // target_res defined
				} // repackable residues
			} else { //redesign/repack interface // inside the interface but outside repackable/redesignable parts
				if ( pose.residue(i).type().is_disulfide_bonded() ) {
					//     task_->nonconst_residue_task(i).prevent_repacking();
					curr_min_sc_[ i ] = false;
					curr_min_bb_[ i ] = false;
				} else if ( automatic_repacking_definition_ ) {
					//     task_->nonconst_residue_task( i ).restrict_to_repacking();
					if ( !min_sc_set() ) curr_min_sc_[ i ] = true;
					if ( !min_bb_set() ) curr_min_bb_[ i ] = true;
				}
			}
		} else { // in interface // non-interface residues are held constant
			if ( automatic_repacking_definition_ ) {
				//    task_->nonconst_residue_task(i).prevent_repacking();
				if ( !min_sc_set() ) curr_min_sc_[ i ] = false;
			}
			if ( !min_bb_set() ) curr_min_bb_[ i ] = true;
		}
	}
	for ( utility::vector1< core::Size >::const_iterator target_it = target_res.begin();
			target_it!=target_res.end(); ++target_it ) {
		// minimize also +-1 amino-acid residues around target residues, unless they are CYD
		if ( !min_sc_set() ) {
			curr_min_sc_[ *target_it ] = true;
			if ( *target_it < pose.size() ) {
				if ( !(pose.residue(*target_it + 1).type().is_disulfide_bonded() ) ) curr_min_sc_[ *target_it + 1 ] = true;
			}
			if ( *target_it > 1 ) {
				if ( !(pose.residue(*target_it - 1).type().is_disulfide_bonded() ) ) curr_min_sc_[ *target_it - 1 ] = true;
			}
		}
		if ( !min_bb_set() ) {
			curr_min_bb_[ *target_it ] = true;
			if ( *target_it < pose.size() ) {
				curr_min_bb_[ *target_it + 1 ] = true;
			}
			if ( *target_it > 1 ) {
				curr_min_bb_[ *target_it - 1 ] = true;
			}
		}
	}
	for ( utility::vector1< core::Size >::const_iterator res=prevent_repacking_.begin(); res!=prevent_repacking_.end(); ++res ) {
		task_->nonconst_residue_task( *res ).prevent_repacking();
	}

	for ( utility::vector1< core::Size >::const_iterator res=restrict_to_repacking_.begin(); res!=restrict_to_repacking_.end(); ++res ) {
		task_->nonconst_residue_task( *res ).restrict_to_repacking();
	}

	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		core::pack::task::parse_resfile(pose, *task_);
	}

	if ( symmetry_ ) {
		core::pack::make_symmetric_PackerTask_by_truncation( in_pose, task_ );
	}
}

utility::tag::XMLSchemaComplexTypeGeneratorOP DesignRepackMover::get_xsd_complex_type() {
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );

	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	rosetta_scripts::attributes_for_parse_score_function( attlist, "scorefxn_repack" );
	rosetta_scripts::attributes_for_parse_score_function( attlist, "scorefxn_minimize" );
	attlist + XMLSchemaAttribute( "repack", xsct_rosetta_bool, "Repack both partners" )
		+ XMLSchemaAttribute( "repack_partner1", xsct_rosetta_bool, "Repack partner 1" )
		+ XMLSchemaAttribute( "repack_partner2", xsct_rosetta_bool, "Repack partner 2" )
		+ XMLSchemaAttribute( "design", xsct_rosetta_bool, "Design both partners" )
		+ XMLSchemaAttribute( "design_partner1", xsct_rosetta_bool, "Design partner 1" )
		+ XMLSchemaAttribute( "design_partner2", xsct_rosetta_bool, "Design partner 2" )
		+ XMLSchemaAttribute::attribute_w_default( "optimize_fold_tree", xsct_rosetta_bool, "Optimize the fold tree", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_rb", xsct_rosetta_bool, "Minimize rigid body degrees of freedom", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_bb", xsct_rosetta_bool, "Minimize backbone of both chains", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_bb_ch1", xsct_rosetta_bool, "Minimize backbone of chain 1", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_bb_ch2", xsct_rosetta_bool, "Minimize backbone of chain 2", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_sc", xsct_rosetta_bool, "Minimize sidechains", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "interface_cutoff_distance", xsct_real, "Interface cutoff distance", "8.0" )
		+ XMLSchemaAttribute::attribute_w_default( "repack_non_ala", xsct_rosetta_bool, "Repack residues that aren't alanine", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "symmetry", xsct_rosetta_bool, "operating on a symmetric pose", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "automatic_repacking_definition", xsct_rosetta_bool, "Is the packing setup automatically?", "1" );

	XMLSchemaSimpleSubelementList ssl;

	AttributeList subtag_attributes;
	subtag_attributes + XMLSchemaAttribute( "pdb_num", xsct_refpose_enabled_residue_number, "Residue number in reference pose or pdb or seqpos numbering" )
		+ XMLSchemaAttribute( "resnum", xsct_refpose_enabled_residue_number, "Residue number in reference pose or pdb or seqpos numbering" );

	ssl.add_simple_subelement( "residue", subtag_attributes, "Tags describing individual residues"/*, 0 minoccurs*/ );

	ct_gen->complex_type_naming_func( & complex_type_name_for_mover );
	ct_gen->add_attributes( attlist );
	ct_gen->add_optional_name_attribute();
	ct_gen->set_subelements_repeatable( ssl );

	return ct_gen;
}

void
DesignRepackMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &data
){
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );

	std::string const scorefxn_repack( protocols::rosetta_scripts::get_score_function_name(tag, "scorefxn_repack" ) );
	std::string const scorefxn_minimize( protocols::rosetta_scripts::get_score_function_name(tag, "scorefxn_minimize" ) );

	if ( !tag->hasOption( "repack_partner1" ) && !tag->hasOption( "repack_partner2" ) && tag->hasOption( "repack" ) ) {
		bool const repack( tag->getOption< bool >( "repack" ) );
		repack_partner1_ = repack_partner2_ = repack;
	} else {
		repack_partner1_ = tag->getOption<bool>( "repack_partner1", true );
		repack_partner2_ = tag->getOption<bool>( "repack_partner2", true );
	}
	if ( !tag->hasOption( "design_partner1" ) && !tag->hasOption( "design_partner2" ) && tag->hasOption( "design" ) ) {
		bool const design( tag->getOption< bool >( "design" ) );
		design_partner1_ = design_partner2_ = design;
	} else {
		design_partner1_ = tag->getOption<bool>( "design_partner1", false );
		design_partner2_ = tag->getOption<bool>( "design_partner2", true );
	}
	TR<<"design partner1: "<<design_partner1_<<'\n';
	TR<<"design partner2: "<<design_partner2_<<'\n';
	TR<<"repack partner1: "<<repack_partner1_<<'\n';
	TR<<"repack partner2: "<<repack_partner2_<<'\n';

	runtime_assert( !( design_partner1_ && !repack_partner1_ ) );
	runtime_assert( !( design_partner2_ && !repack_partner2_ ) );
	optimize_foldtree_ = tag->getOption<bool>( "optimize_fold_tree", true );
	if ( tag->hasOption( "minimize_rb" ) ) {
		min_rb( tag->getOption<bool>( "minimize_rb", true ) );
	}

	if ( tag->hasOption( "minimize_bb" ) ) {
		if ( tag->getOption<bool>( "minimize_bb", true ) ) {
			min_bb( utility::pointer::make_shared< core::select::residue_selector::TrueResidueSelector >() );
		} else {
			min_bb( utility::pointer::make_shared< core::select::residue_selector::FalseResidueSelector >() );
		}
	}

	if ( tag->hasOption( "minimize_bb_ch1" ) ||  ( tag->hasOption( "minimize_bb_ch2" ) ) ) {
		if ( tag->hasOption( "minimize_bb" ) ) {
			TR.Warning << "Setting minimize_bb_ch1 or minimize_bb_ch2 for DesignRepackMover causes it to ignore minimize_bb settings." << std::endl;
		}
		core::select::residue_selector::ResidueSelectorOP ch1;
		core::select::residue_selector::ResidueSelectorOP ch2;
		if ( tag->getOption<bool>( "minimize_bb_ch1", false ) ) {
			ch1 = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >(1);
		}
		if ( tag->getOption<bool>( "minimize_bb_ch2", false ) ) {
			ch2 = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >(2);
		}
		if ( ch1 != nullptr && ch2 != nullptr ) {
			min_bb( core::select::residue_selector::OR_combine( ch1, ch2 ) );
		} else if ( ch1 != nullptr ) {
			min_bb( ch1 );
		} else if ( ch2 != nullptr ) {
			min_bb( ch2 );
		} else {
			min_bb( utility::pointer::make_shared< core::select::residue_selector::FalseResidueSelector >() );
		}
	}//end specificatino of bb mininization

	if ( tag->hasOption( "minimize_sc" ) ) {
		if ( tag->getOption<bool>( "minimize_sc", true ) ) {
			min_sc( utility::pointer::make_shared< core::select::residue_selector::TrueResidueSelector >() );
		} else {
			min_sc( utility::pointer::make_shared< core::select::residue_selector::FalseResidueSelector >() );
		}
	}
	interface_distance_cutoff_ = tag->getOption<core::Real>( "interface_cutoff_distance", 8.0 );
	utility::vector0< TagCOP > const & repack_tags( tag->getTags() );

	utility::vector1< std::string > target_res_vec;
	for ( auto repack_ptr : repack_tags ) {
		if ( repack_ptr->getName() == "residue" ) {
			target_res_vec.push_back( core::pose::get_resnum_string( repack_ptr ) );
		}
	}
	if ( target_res_vec.size() ) {
		target_residues( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( utility::join(target_res_vec,",") ) );
	}

	repack_non_ala_ = tag->getOption<bool>( "repack_non_ala", true );

	symmetry_ = tag->getOption< bool >( "symmetry", false );

	using namespace core::scoring;
	scorefxn_repack_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_repack", data )->clone();
	scorefxn_minimize_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_minimize", data )->clone();
	automatic_repacking_definition_ = tag->getOption<bool>( "automatic_repacking_definition", true );
	TR<<"repack scorefxn "<<scorefxn_repack<<" and minimize scorefxn "<<scorefxn_minimize<<" automatic_repacking_definition set to "<<automatic_repacking_definition_<<" optimize fold tree="<<optimize_foldtree_<< std::endl;
}

void
DesignRepackMover::set_scorefxn_repack( core::scoring::ScoreFunctionCOP scorefxn ) {
	scorefxn_repack_ = scorefxn->clone();
}

void
DesignRepackMover::set_scorefxn_minimize( core::scoring::ScoreFunctionCOP scorefxn ) {
	scorefxn_minimize_ = scorefxn->clone();
}

void
DesignRepackMover::clear_task_factory(){
	if ( task_factory_ ) {
		task_factory_->clear();
	}
}

core::pack::task::PackerTaskCOP
DesignRepackMover::task() const {
	return task_;
}

core::pack::task::PackerTaskOP &
DesignRepackMover::task() {
	return task_;
}

core::pack::task::TaskFactoryOP &
DesignRepackMover::task_factory() {
	return task_factory_;
}

core::pack::task::TaskFactoryOP
DesignRepackMover::task_factory() const {
	return task_factory_;
}

DesignRepackMover::~DesignRepackMover() = default;

core::scoring::ScoreFunctionOP
DesignRepackMover::scorefxn_repack() const {
	return scorefxn_repack_;
}

core::scoring::ScoreFunctionOP
DesignRepackMover::scorefxn_minimize() const {
	return scorefxn_minimize_;
}

void
DesignRepackMover::target_residues( core::select::residue_selector::ResidueSelectorCOP setting ) {
	target_residues_ = setting;
}

core::select::residue_selector::ResidueSelectorCOP
DesignRepackMover::target_residues() const {
	return target_residues_;
}

utility::vector1< core::Size >
DesignRepackMover::target_residues( core::pose::Pose const & pose ) const {
	utility::vector1< core::Size > target_res;
	if ( target_residues_ ) {
		target_res = core::select::get_residues_from_subset( target_residues_->apply( pose ) );
		if ( target_res.empty() ) {
			TR.Warning << "A residue selection was given, but it's empty!" << std::endl;
		}
	}
	return target_res;
}

} //simple_moves
} //protocols
