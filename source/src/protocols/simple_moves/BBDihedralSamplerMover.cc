// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/BBDihedralSamplerMover.cc
/// @brief Mover interface to BBDihedralSampler.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/simple_moves/BBDihedralSamplerMover.hh>
#include <protocols/simple_moves/BBDihedralSamplerMoverCreator.hh>
#include <protocols/simple_moves/bb_sampler/BBDihedralSampler.hh>

#include <core/pose/Pose.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
#include <core/select/residue_selector/util.hh>


#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <numeric/random/random.hh>

#include <utility/excn/Exceptions.hh> // AUTO IWYU For Exception
#include <utility/stream_util.hh> // AUTO IWYU For operator<<

static basic::Tracer TR( "protocols.simple_moves.BBDihedralSamplerMover" );


namespace protocols {
namespace simple_moves {
using namespace protocols::simple_moves::bb_sampler;
using namespace core::kinematics;

BBDihedralSamplerMover::BBDihedralSamplerMover():
	protocols::moves::Mover( "BBDihedralSamplerMover" )
{

}

BBDihedralSamplerMover::BBDihedralSamplerMover( BBDihedralSamplerOP sampler):
	protocols::moves::Mover( "BBDihedralSamplerMover" )

{
	set_sampler( sampler );
}

BBDihedralSamplerMover::BBDihedralSamplerMover( BBDihedralSamplerOP sampler, core::select::residue_selector::ResidueSelectorCOP selector):
	protocols::moves::Mover( "BBDihedralSamplerMover" )
{
	set_sampler( sampler );
	set_residue_selector( selector );
}

BBDihedralSamplerMover::~BBDihedralSamplerMover()= default;

BBDihedralSamplerMover::BBDihedralSamplerMover( BBDihedralSamplerMover const & src ):
	protocols::moves::Mover( src ),
	sampler_torsion_types_(src.sampler_torsion_types_),
	bb_residues_(src.bb_residues_),
	sampler_torsions_(src.sampler_torsions_)
{
	samplers_.clear();
	for ( auto const & kv : src.samplers_ ) {
		utility::vector1< BBDihedralSamplerOP > sampler_list;
		for ( BBDihedralSamplerOP sampler : kv.second ) {
			BBDihedralSamplerOP new_sampler = sampler->clone();
			sampler_list.push_back( new_sampler );
		}
		samplers_[kv.first] = sampler_list;
	}
	if ( src.selector_ ) selector_ = src.selector_->clone();
}

void
BBDihedralSamplerMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
)
{

}

protocols::moves::MoverOP
BBDihedralSamplerMover::clone() const{
	return utility::pointer::make_shared< BBDihedralSamplerMover >( *this );
}


moves::MoverOP
BBDihedralSamplerMover::fresh_instance() const
{
	return utility::pointer::make_shared< BBDihedralSamplerMover >();
}

std::string
BBDihedralSamplerMover::get_name() const {
	return "BBDihedralSamplerMover";
}

void
BBDihedralSamplerMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

std::ostream &operator<< (std::ostream &os, BBDihedralSamplerMover const &mover)
{
	mover.show(os);
	return os;
}

void
BBDihedralSamplerMover::set_name_of_last_sampler_used( std::string name_of_last_sampler_used ){
	name_of_last_sampler_used_ = name_of_last_sampler_used;
}

std::string
BBDihedralSamplerMover::get_name_of_last_sampler_used(){
	return name_of_last_sampler_used_;
}


void
BBDihedralSamplerMover::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector->clone();
	sampler_torsions_.clear();
	bb_residues_.clear();

}

void
BBDihedralSamplerMover::set_single_resnum( core::Size resnum ){
	bb_residues_.clear();
	sampler_torsions_.clear();
	bb_residues_.push_back( resnum );
}

void
BBDihedralSamplerMover::set_sampler( bb_sampler::BBDihedralSamplerOP sampler ){
	samplers_.clear();
	sampler_torsion_types_.clear();
	sampler_torsions_.clear();
	samplers_[ sampler->get_torsion_type() ]; //Initialize vector
	samplers_[ sampler->get_torsion_type() ].push_back( sampler );


	sampler_torsion_types_.push_back( sampler->get_torsion_type() );

}

void
BBDihedralSamplerMover::add_sampler( bb_sampler::BBDihedralSamplerOP sampler ){
	sampler_torsions_.clear();
	if ( samplers_.count( sampler->get_torsion_type() ) != 0 ) {
		samplers_[ sampler->get_torsion_type() ].push_back( sampler );
		//sampler_torsion_types_.push_back( sampler->get_torsion_type() );
	} else {
		samplers_[ sampler->get_torsion_type() ];
		samplers_[ sampler->get_torsion_type() ].push_back( sampler );

		sampler_torsion_types_.push_back( sampler->get_torsion_type() );
	}

}

void
BBDihedralSamplerMover::setup_all_bb_residues( core::pose::Pose const & pose) {

	using namespace core::select::residue_selector;
	if ( TR.Debug.visible() ) {
		TR.Debug << "Setting up all bb residues" << std::endl;
	}

	utility::vector1< bool > subset( pose.total_residue(), true);
	selector_ = utility::pointer::make_shared< ReturnResidueSubsetSelector >( subset );
	bb_residues_ = core::select::residue_selector::selection_positions(subset);
}

///@brief Limit the torsions that are sampled using a mask.
/// Mask is a vector of torsions we are allowed to sample at each residue.
///  NOT every residue must be present - only the ones we need to mask.
///  This helps for carbohydrates, since not every residue has all possible dihedrals, even if we have samplers to use.
void
BBDihedralSamplerMover::set_dihedral_mask( std::map< core::Size, utility::vector1< core::Size >> mask ){
	sampler_torsions_.clear();
	dihedral_mask_ = mask;
	TR.Debug << "Using dihedral mask: " << mask << std::endl;
}

void
BBDihedralSamplerMover::setup_samplers( core::pose::Pose const &  ) {

	if ( TR.Debug.visible() ) {
		TR.Debug << "Initializing samplers" << std::endl;
		TR.Debug << "Overall torsion types available: " << sampler_torsion_types_ << std::endl;
	}

	//TR.Info << "Torsion Types: " << utility::to_string( sampler_torsion_types_) << std::endl;

	//Apply any mask set for the residue, otherwise set the torsion from the samplers we have.
	// For every residue specified by either 1) set_single_resnum, 2) set_residue_selector, or 3) setup_all_bb_residues
	// i.e. for all residues set to be available for sampling
	for ( core::Size resnum : bb_residues_ ) {
		// The resnum to vector of torsion IDs mapping should at first be empty
		sampler_torsions_[ resnum ];
		//TR.Info << "Resnum " << resnum << std::endl;
		// For every dihedral available  based on the torsion types specified by all given samplers...
		// (Note that every residue may not have the same torsion types available to it)
		// Example: some sugar residues have only phi and psi while another may also have omega
		for ( core::Size dihedral : sampler_torsion_types_ ) {
			// If the dihedral_mask_ contains this residue (note not every residue needs a mask)
			if ( dihedral_mask_.count( resnum ) ) {
				// And if the dihedral_mask_ allows the dihedral for this residue
				if ( dihedral_mask_[resnum].contains( dihedral ) ) {
					// Then keep it as an available dihedral to sample for this residue
					sampler_torsions_[ resnum ].push_back( dihedral );
				}
			} else {
				// If the dihedral_mask_ does not contain this residue, then assume this is an available torsion for sampling
				sampler_torsions_[ resnum ].push_back( dihedral );
			}
		}
	}
	// Mapping of residue number to a vector of torsion IDs to sample on
	// Ex. [ [100] : [1, 2], [101] : [1, 2, 3], ... ]
	TR.Debug << "Sampling torsions available: " << sampler_torsions_ << std::endl;
}

void
BBDihedralSamplerMover::apply( core::pose::Pose & pose ){
	utility::vector1< bool > subset;
	if ( selector_ && bb_residues_.size() == 0 ) {
		TR.Debug << "Setting up from selector" << std::endl;
		subset = selector_->apply(pose);
		bb_residues_ = core::select::residue_selector::selection_positions(subset);

		if ( bb_residues_.size() == 0 ) {
			TR << "No BB residues to model (remember - no data for root!) Returning) " << std::endl;
			set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
			return;
		}
	}

	if ( samplers_.size() == 0 ) {
		utility_exit_with_message(" No Sampler set for BBDihedralSamplerMover!");
	}

	if ( bb_residues_.size() == 0 ) {
		setup_all_bb_residues( pose );
	}

	if ( sampler_torsions_.empty() ) {
		setup_samplers( pose );
	}

	if ( TR.Debug.visible() ) { TR.Debug << "Available residues: " <<
		utility::to_string( bb_residues_ ) << std::endl;
	}

	//TR.Info << "Selecting a residue number" << std::endl;
	core::Size const resnum_index =
		numeric::random::rg().random_range2( 1, bb_residues_.size() );
	core::Size const resnum = bb_residues_[ resnum_index ];
	if ( TR.Debug.visible() ) { TR.Debug << "Selected residue " <<
		resnum << " at index " << resnum_index << std::endl;
	}

	//Get the one sampler or choose a sampler:
	bb_sampler::BBDihedralSamplerCOP sampler;
	if ( sampler_torsions_[ resnum ].size() == 0 ) {
		utility_exit_with_message(" BBDihedralSamplerMover - a chosen "
			"residue number has no dihedral union "
			"between movemaps and set bb samplers. "
			"We should never be here!");
	}

	//Choose a TorsionType based on what kind of samplers we have and the MoveMap of the residue.
	core::Size const torsion_index =
		numeric::random::rg().random_range( 1, sampler_torsions_[ resnum ].size() );
	//TR.Info << "Torsion Index " << torsion_index << std::endl;
	core::Size const torsion = sampler_torsions_[ resnum ][ torsion_index ];

	//TR.Info << "Selected torsion " << torsion << std::endl;
	//Choose a sampler from the samplers available for that torsion type.
	core::Size const sampler_index =
		numeric::random::rg().random_range( 1, samplers_[ torsion ].size() );
	sampler = samplers_[ torsion ][ sampler_index ];

	// Apply the sampler
	TR.Debug << "Optimizing "<< resnum << " with " <<
		sampler->get_name() << " at torsion " << torsion << std::endl;
	// And store the name of the sampler as the one last used
	set_name_of_last_sampler_used( sampler->get_name() );

	try {
		sampler->set_torsion_to_pose( pose, resnum );
		set_last_move_status(protocols::moves::MS_SUCCESS);

	} catch ( utility::excn::Exception& excn ) {
		TR.Error << "Could not set torsion for resnum " << resnum << std::endl;
		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
	}

}


/////////////// Creator ///////////////

protocols::moves::MoverOP
BBDihedralSamplerMoverCreator::create_mover() const {
	return utility::pointer::make_shared< BBDihedralSamplerMover >();
}

std::string
BBDihedralSamplerMoverCreator::keyname() const {
	return BBDihedralSamplerMoverCreator::mover_name();
}

std::string
BBDihedralSamplerMoverCreator::mover_name(){
	return "BBDihedralSamplerMover";
}

/// @brief Provide the citation.
void
BBDihedralSamplerMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	using namespace basic::citation_manager;

	citations.add(
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		get_name(),
		CitedModuleType::Mover,
		"Jared Adolf-Bryfogle",
		"The Scripps Research Institute, La Jolla, CA",
		"jadolfbr@gmail.com"
		)
	);

	citations.add( selector_ );
}

} //protocols
} //carbohydrates


