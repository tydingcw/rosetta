// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/RNA_DataBackboneEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das


// Unit headers
#include <core/energy_methods/RNA_DataBackboneEnergy.hh>
#include <core/energy_methods/RNA_DataBackboneEnergyCreator.hh>

// Package headers
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/pose/rna/RNA_DataInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/id/AtomID.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// Utility headers
#include <ObjexxFCL/FArray1D.hh>


// C++
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>

using namespace core::chemical::rna;

namespace core {
namespace energy_methods {

using Vector = numeric::xyzVector<core::Real>;

/// @details This must return a fresh instance of the RNA_DataBackboneEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
RNA_DataBackboneEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< RNA_DataBackboneEnergy >();
}

core::scoring::ScoreTypes
RNA_DataBackboneEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( rna_data_backbone );
	return sts;
}


/// c-tor
RNA_DataBackboneEnergy::RNA_DataBackboneEnergy() :
	parent( utility::pointer::make_shared< RNA_DataBackboneEnergyCreator >() ),
	dist_cutoff_( 9.0 ),
	dist_fade_( 1.0  ),
	well_depth_burial_( -0.05 ),
	well_depth_exposed_( 0.01 ),
	burial_function_( utility::pointer::make_shared< core::scoring::func::FadeFunc >( -10.0 /*cutoff_lower*/, dist_cutoff_, dist_fade_, 1.0 /*well_depth*/ ) )
{
	initialize_atom_numbers_sugar();
}


//clone
core::scoring::methods::EnergyMethodOP
RNA_DataBackboneEnergy::clone() const
{
	return utility::pointer::make_shared< RNA_DataBackboneEnergy >();
}

///////////////////////////////////////////////////////////
void
RNA_DataBackboneEnergy::initialize_atom_numbers_sugar() {

	using namespace core::chemical;
	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	ResidueTypeCOP const & rsd_type( rsd_set->get_representative_type_aa( na_rad ) ); //Check out adenine.

	atom_numbers_sugar_.clear();
	atom_numbers_sugar_.push_back( rsd_type->atom_index( "C1'" ) );
	atom_numbers_sugar_.push_back( rsd_type->atom_index( "C2'" ) );
	atom_numbers_sugar_.push_back( rsd_type->atom_index( "C3'" ) );
	atom_numbers_sugar_.push_back( rsd_type->atom_index( "C4'" ) );
	atom_numbers_sugar_.push_back( rsd_type->atom_index( "C5'" ) );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
RNA_DataBackboneEnergy::setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	// need to make sure scoring info is in there...
	core::scoring::rna::nonconst_rna_scoring_info_from_pose( pose );

	pose.update_residue_neighbors();
}


void
RNA_DataBackboneEnergy::setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}

/////////////////////////////////////////////////////////////////
void
RNA_DataBackboneEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	pose.update_residue_neighbors();
}


//////////////////////////////////////////////////////////////////////////////////////////
//
//
//
void
RNA_DataBackboneEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const {

	if ( !rsd1.is_RNA() ) return;
	if ( !rsd2.is_RNA() ) return;
	if ( rsd1.is_virtual_residue() ) return;
	if ( rsd2.is_virtual_residue() ) return;
	if ( rsd1.has_variant_type( chemical::REPLONLY ) ) return;
	if ( rsd2.has_variant_type( chemical::REPLONLY ) ) return;

	//  rna_filtered_base_base_info.set_calculated( false );
	core::scoring::rna::RNA_ScoringInfo const & rna_scoring_info( core::scoring::rna::rna_scoring_info_from_pose( pose ) );
	pose::rna::RNA_DataInfo const & rna_data_info( rna_scoring_info.rna_data_info() );
	ObjexxFCL::FArray1D < bool > const & rna_data_backbone_burial( rna_data_info.backbone_burial() );
	ObjexxFCL::FArray1D < bool > const & rna_data_backbone_exposed( rna_data_info.backbone_exposed() );

	if ( rna_data_backbone_burial.size() == 0.0 ) return;

	debug_assert( rna_data_backbone_burial.size() == pose.size() );

	if ( rna_data_backbone_burial( rsd1.seqpos() ) ) {
		emap[ core::scoring::rna_data_backbone ]         += well_depth_burial_ * get_sugar_env_score( rsd1 /*buried sugar*/, rsd2 /*other*/ );
	} else if ( rna_data_backbone_exposed( rsd1.seqpos() )  ) {
		emap[ core::scoring::rna_data_backbone ]         += well_depth_exposed_ * get_sugar_env_score( rsd1 /*buried sugar*/, rsd2 /*other*/ );
	}


	if ( rna_data_backbone_burial( rsd2.seqpos() ) ) {
		emap[ core::scoring::rna_data_backbone ]         += well_depth_burial_ * get_sugar_env_score( rsd2 /*buried sugar*/, rsd1 /*other*/ );
	} else if ( rna_data_backbone_exposed( rsd2.seqpos() )  ) {
		emap[ core::scoring::rna_data_backbone ]         += well_depth_exposed_ * get_sugar_env_score( rsd2 /*buried sugar*/, rsd1 /*other*/ );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
Real
RNA_DataBackboneEnergy::get_sugar_env_score( core::conformation::Residue const & rsd_buried, core::conformation::Residue const & rsd_other ) const
{
	Real score( 0.0 );
	bool is_coarse = ( rsd_buried.is_coarse() );

	utility::vector1< Size > const & atom_numbers_sugar = ( is_coarse ? atom_numbers_sugar_coarse_ : atom_numbers_sugar_ );
	Size const num_sugar_atoms( atom_numbers_sugar.size() );

	for ( Size j = 1; j <= num_sugar_atoms; j++ ) {
		numeric::xyzVector< Real > const & sugar_atom( rsd_buried.xyz( atom_numbers_sugar[j] ) );
		for ( Size m = 1; m <= rsd_other.nheavyatoms(); m++ ) {
			numeric::xyzVector< Real > const & other_atom( rsd_other.xyz( m ) );
			Real const dist = ( other_atom - sugar_atom ).length();
			Real const burial_score = burial_function_->func( dist );
			score += burial_score;
		}
	}

	return score;
}


// Need to define atom derivative!!
///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
bool
RNA_DataBackboneEnergy::check_sugar_atom( Size const & n ) const {

	for ( Size i = 1; i <= atom_numbers_sugar_.size(); i++ ) {
		if ( atom_numbers_sugar_[i] == n ) return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////
void
RNA_DataBackboneEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( i ) );

	if ( rsd1.is_virtual_residue() ) return;
	if ( rsd1.has_variant_type( chemical::REPLONLY ) ) return;
	if ( rsd1.is_repulsive( m ) ) return;

	if ( ! rsd1.is_RNA() ) return;
	if ( m > rsd1.nheavyatoms() ) return;

	core::scoring::rna::RNA_ScoringInfo  const & rna_scoring_info( core::scoring::rna::rna_scoring_info_from_pose( pose ) );
	pose::rna::RNA_DataInfo const & rna_data_info( rna_scoring_info.rna_data_info() );
	ObjexxFCL::FArray1D < bool > const & rna_data_backbone_burial( rna_data_info.backbone_burial() );
	ObjexxFCL::FArray1D < bool > const & rna_data_backbone_exposed( rna_data_info.backbone_exposed() );

	if ( rna_data_backbone_burial.size() == 0.0 ) return;

	debug_assert( rna_data_backbone_burial.size() == pose.size() );

	Vector const heavy_atom_i( rsd1.xyz( m ) );
	bool const pos1_fixed( domain_map( i ) != 0 );

	// cached energies object
	core::scoring::Energies const & energies( pose.energies() );

	// the neighbor/energy links
	core::scoring::EnergyGraph const & energy_graph( energies.energy_graph() );

	for ( utility::graph::Graph::EdgeListConstIter
			iter  = energy_graph.get_node( i )->const_edge_list_begin(),
			itere = energy_graph.get_node( i )->const_edge_list_end();
			iter != itere; ++iter ) {

		Size const j( ( *iter )->get_other_ind( i ) );
		if ( pos1_fixed && domain_map( i ) == domain_map( j ) ) continue; //Fixed w.r.t. one another.

		conformation::Residue const & rsd2( pose.residue( j ) );
		if ( ! rsd2.is_RNA() ) continue;

		if ( rsd2.is_virtual_residue() ) continue;
		if ( rsd2.has_variant_type( chemical::REPLONLY ) ) continue;

		// This could be faster if split into separate loops.
		if ( rna_data_backbone_burial( i ) && check_sugar_atom( m ) ) { // other heavy atoms are possible buriers
			for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {
				if ( rsd2.is_repulsive( n ) ) continue;

				Vector const heavy_atom_j( rsd2.xyz( n ) );
				Vector r = heavy_atom_j - heavy_atom_i;
				Real const dist = r.length();
				Real const deriv = burial_function_->dfunc( dist );

				Vector const force_vector_i = deriv * r / dist;
				Vector const f1 = cross( force_vector_i, heavy_atom_j );
				Vector const f2 = force_vector_i;
				F1 +=  -1.0 * well_depth_burial_ * weights[ core::scoring::rna_data_backbone ] * f1;
				F2 +=  -1.0 * well_depth_burial_ * weights[ core::scoring::rna_data_backbone ] * f2;
			}
		}

		if ( rna_data_backbone_exposed( i ) && check_sugar_atom( m ) ) { // other heavy atoms are possible buriers
			for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {
				if ( rsd2.is_repulsive( n ) ) continue;

				Vector const heavy_atom_j( rsd2.xyz( n ) );
				Vector r = heavy_atom_j - heavy_atom_i;
				Real const dist = r.length();
				Real const deriv = burial_function_->dfunc( dist );

				Vector const force_vector_i = deriv * r / dist;
				Vector const f1 = cross( force_vector_i, heavy_atom_j );
				Vector const f2 = force_vector_i;
				F1 +=  -1.0 * well_depth_exposed_ * weights[ core::scoring::rna_data_backbone ] * f1;
				F2 +=  -1.0 * well_depth_exposed_ * weights[ core::scoring::rna_data_backbone ] * f2;
			}
		}


		if ( rna_data_backbone_burial( j ) ) { // other residue's sugar atoms might be buried.
			for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {
				if ( ! check_sugar_atom( n ) ) continue;
				if ( rsd2.is_repulsive( n ) ) continue;

				Vector const heavy_atom_j( rsd2.xyz( n ) );
				Vector r = heavy_atom_j - heavy_atom_i;
				Real const dist = r.length();
				Real const deriv = burial_function_->dfunc( dist );

				Vector const force_vector_j = deriv * (  - r ) / dist;
				Vector const f1 = -1.0 * cross( force_vector_j, heavy_atom_i );
				Vector const f2 = -1.0 * force_vector_j;
				F1 +=  -1.0 * well_depth_burial_ * weights[ core::scoring::rna_data_backbone ] * f1;
				F2 +=  -1.0 * well_depth_burial_ * weights[ core::scoring::rna_data_backbone ] * f2;
			}
		}


		if ( rna_data_backbone_exposed( j ) ) { // other residue's sugar atoms might be buried.
			for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {
				if ( ! check_sugar_atom( n ) ) continue;
				if ( rsd2.is_repulsive( n ) ) continue;

				Vector const heavy_atom_j( rsd2.xyz( n ) );
				Vector r = heavy_atom_j - heavy_atom_i;
				Real const dist = r.length();
				Real const deriv = burial_function_->dfunc( dist );

				Vector const force_vector_j = deriv * (  - r ) / dist;
				Vector const f1 = -1.0 * cross( force_vector_j, heavy_atom_i );
				Vector const f2 = -1.0 * force_vector_j;
				F1 +=  -1.0 * well_depth_exposed_ * weights[ core::scoring::rna_data_backbone ] * f1;
				F2 +=  -1.0 * well_depth_exposed_ * weights[ core::scoring::rna_data_backbone ] * f2;
			}
		}
	}
}


/// @brief RNA_DataBackboneEnergy distance cutoff
Distance
RNA_DataBackboneEnergy::atomic_interaction_cutoff() const
{
	return 0.0; /// Uh, I don't know.
}
core::Size
RNA_DataBackboneEnergy::version() const
{
	return 1; // Initial versioning
}

} //scoring
} //core
