// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/disulfides/FullatomDisulfideEnergy.cc
/// @brief  Disulfide Energy class implementation
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/energy_methods/FullatomDisulfideEnergy.hh>
#include <core/energy_methods/FullatomDisulfideEnergyCreator.hh>

// Package headers
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.hh>
#include <core/scoring/disulfides/FullatomDisulfideEnergyContainer.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/methods/Methods.hh>
#include <basic/datacache/CacheableData.hh>

#include <core/scoring/disulfides/DisulfideAtomIndices.hh>
#include <utility/vector1.hh>


namespace core {
namespace energy_methods {

class DisulfMinData : public basic::datacache::CacheableData
{
public:
	using CacheableDataOP = basic::datacache::CacheableDataOP;
public:
	DisulfMinData( conformation::Residue const & res1, conformation::Residue const & res2 );
	~DisulfMinData() override;
	CacheableDataOP clone() const override;

	/// @brief which_res should be 1 or 2
	void set_res_inds( Size which_res, core::scoring::disulfides::DisulfideAtomIndices const & dais );
	core::scoring::disulfides::DisulfideAtomIndices const & res_inds( Size which_res ) const;

private:
	core::scoring::disulfides::DisulfideAtomIndices res1_inds_;
	core::scoring::disulfides::DisulfideAtomIndices res2_inds_;
};

using DisulfMinDataOP = utility::pointer::shared_ptr<DisulfMinData>;
using DisulfMinDataCOP = utility::pointer::shared_ptr<const DisulfMinData>;


DisulfMinData::DisulfMinData(
	conformation::Residue const & res1,
	conformation::Residue const & res2
) :
	res1_inds_( res1 ),
	res2_inds_( res2 )
{}

DisulfMinData::~DisulfMinData() = default;

DisulfMinData::CacheableDataOP DisulfMinData::clone() const { return utility::pointer::make_shared< DisulfMinData >( *this ); }

/// @brief which_res should be 1 or 2
void DisulfMinData::set_res_inds( Size which_res, core::scoring::disulfides::DisulfideAtomIndices const & dais )
{
	debug_assert( which_res == 1 || which_res == 2 );
	if ( which_res == 1 ) res1_inds_ = dais;
	else res2_inds_ = dais;
}

core::scoring::disulfides::DisulfideAtomIndices const &
DisulfMinData::res_inds( Size which_res ) const
{
	debug_assert( which_res == 1 || which_res == 2 );
	return which_res == 1 ? res1_inds_ : res2_inds_;
}


/// @details This must return a fresh instance of the FullatomDisulfideEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
FullatomDisulfideEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< FullatomDisulfideEnergy >( core::scoring::ScoringManager::get_instance()->get_FullatomDisulfidePotential() );
}

core::scoring::ScoreTypes
FullatomDisulfideEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( dslf_ss_dst );
	sts.push_back( dslf_cs_ang );
	sts.push_back( dslf_ss_dih );
	sts.push_back( dslf_ca_dih );
	sts.push_back( dslf_cbs_ds );
	sts.push_back( dslf_fa13 );
	return sts;
}


FullatomDisulfideEnergy::FullatomDisulfideEnergy( core::scoring::disulfides::FullatomDisulfidePotential const & potential )
:
	parent( utility::pointer::make_shared< FullatomDisulfideEnergyCreator >() ),
	potential_( potential )
{}

FullatomDisulfideEnergy::~FullatomDisulfideEnergy() = default;

// EnergyMethod Methods:

core::scoring::methods::EnergyMethodOP
FullatomDisulfideEnergy::clone() const
{
	return utility::pointer::make_shared< FullatomDisulfideEnergy >( potential_ );
}

void
FullatomDisulfideEnergy::ensure_lrenergy_container_is_up_to_date(
	pose::Pose & pose
) const
{
	using namespace core::scoring::methods;
	using namespace core::scoring::disulfides;

	if ( pose.energies().long_range_container( fa_disulfide_energy ) == nullptr ) {
		FullatomDisulfideEnergyContainerOP dec = utility::pointer::make_shared< FullatomDisulfideEnergyContainer >( pose );
		pose.energies().set_long_range_container( fa_disulfide_energy, dec );
	} else {
		FullatomDisulfideEnergyContainerOP dec = FullatomDisulfideEnergyContainerOP (
			utility::pointer::static_pointer_cast< FullatomDisulfideEnergyContainer > ( pose.energies().nonconst_long_range_container( fa_disulfide_energy ) ) );
		dec->update( pose );
		if ( dec->num_residues() != pose.conformation().size() ) {
			FullatomDisulfideEnergyContainerOP dec2 = utility::pointer::make_shared< FullatomDisulfideEnergyContainer >( pose );
			pose.energies().set_long_range_container( fa_disulfide_energy, dec2 );
		}
	}
}

void
FullatomDisulfideEnergy::setup_for_scoring(
	pose::Pose & pose,
	core::scoring::ScoreFunction const & ) const
{
	ensure_lrenergy_container_is_up_to_date( pose );
}

void
FullatomDisulfideEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const & ,
	utility::vector1< bool > const &
) const
{
	ensure_lrenergy_container_is_up_to_date( pose );
}

/// @details returns true if both residues are cys, if both are disulfide-cys, and then
/// if all of these conditions have been satisfied, if residue1's SG atom connects to residue 2,
/// and if residue 2's SG atom connects to residue 1
/// @note Correction to above: returns true if both residues are possible disulfide-forming
/// residues, and both share a disulfide bond.
bool
FullatomDisulfideEnergy::defines_score_for_residue_pair(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	bool res_moving_wrt_eachother
) const
{
	using namespace chemical;
	std::string atom1 = res1.type().get_disulfide_atom_name();
	std::string atom2 = res2.type().get_disulfide_atom_name();

	// Formerly tested if the residue had its own disulfide atom name,
	// but checking if it is NONE is just clearer (and heads off the
	// case where a new residue is added with an atom named NONE
	return res_moving_wrt_eachother
		&& res1.has_variant_type( DISULFIDE ) && res2.has_variant_type( DISULFIDE )
		&& res1.type().get_disulfide_atom_name() != "NONE"
		&& res1.connect_map( res1.type().residue_connection_id_for_atom( res1.atom_index( atom1 ) ) ).resid() == res2.seqpos()
		&& res2.type().get_disulfide_atom_name() != "NONE"
		&& res2.connect_map( res2.type().residue_connection_id_for_atom( res2.atom_index( atom2 ) ) ).resid() == res1.seqpos();
}

bool
FullatomDisulfideEnergy::use_extended_residue_pair_energy_interface() const { return true; }

bool
FullatomDisulfideEnergy::minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

void
FullatomDisulfideEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::scoring::ResPairMinimizationData const & min_data,
	pose::Pose const &,
	core::scoring::ScoreFunction const &sfxn,
	core::scoring::EnergyMap & emap
) const
{
	using namespace core::scoring::disulfides;

	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	// Ignore VIRTUAL residues.
	if ( rsd1.is_virtual_residue() || rsd2.is_virtual_residue() ) {
		return;
	}

	conformation::Residue const & rsdl( rsd1.seqpos() < rsd2.seqpos() ? rsd1 : rsd2 );
	conformation::Residue const & rsdu( rsd1.seqpos() < rsd2.seqpos() ? rsd2 : rsd1 );

	debug_assert( utility::pointer::dynamic_pointer_cast< DisulfMinData const > ( min_data.get_data( core::scoring::fa_dslf_respair_data ) ) );
	DisulfMinData const & disulf_inds = static_cast< DisulfMinData const & > ( *min_data.get_data( core::scoring::fa_dslf_respair_data ) );

	//fpd old version
	if ( sfxn.has_nonzero_weight(core::scoring::dslf_ss_dst) || sfxn.has_nonzero_weight(core::scoring::dslf_cs_ang) || sfxn.has_nonzero_weight(core::scoring::dslf_ss_dih) || sfxn.has_nonzero_weight(core::scoring::dslf_ca_dih) ) {
		Energy distance_score_this_disulfide;
		Energy csangles_score_this_disulfide;
		Energy dihedral_score_this_disulfide;
		Energy ca_dihedral_sc_this_disulf;
		bool truefalse_fa_disulf;

		potential_.score_this_disulfide_old(
			rsdl, rsdu,
			disulf_inds.res_inds( 1 ),
			disulf_inds.res_inds( 2 ),
			distance_score_this_disulfide,
			csangles_score_this_disulfide,
			dihedral_score_this_disulfide,
			ca_dihedral_sc_this_disulf,
			truefalse_fa_disulf
		);

		//amw
		/*std::cout << "scoring dslf old: " <<std::endl;
		std::cout << emap[ core::scoring::dslf_ss_dst ] << " += " << distance_score_this_disulfide << std::endl;
		std::cout << emap[ core::scoring::dslf_cs_ang ] << " += " << csangles_score_this_disulfide << std::endl;
		std::cout << emap[ core::scoring::dslf_ss_dih ] << " += " << dihedral_score_this_disulfide << std::endl;
		std::cout << emap[ core::scoring::dslf_ca_dih ] << " += " << ca_dihedral_sc_this_disulf << std::endl;
		*/

		emap[ core::scoring::dslf_ss_dst ] += distance_score_this_disulfide;
		emap[ core::scoring::dslf_cs_ang ] += csangles_score_this_disulfide;
		emap[ core::scoring::dslf_ss_dih ] += dihedral_score_this_disulfide;
		emap[ core::scoring::dslf_ca_dih ] += ca_dihedral_sc_this_disulf;
	}

	//fpd new version
	if ( sfxn.has_nonzero_weight( core::scoring::dslf_fa13 ) ) {
		Energy score_i;
		potential_.score_this_disulfide(
			rsdl, rsdu,
			disulf_inds.res_inds( 1 ), disulf_inds.res_inds( 2 ),
			score_i );
		//amw
		/*std::cout << "scoring dslf new: " <<std::endl;
		std::cout << emap[ core::scoring::dslf_fa13 ] << " += " << score_i << std::endl;
		*/
		emap[ core::scoring::dslf_fa13 ] += score_i;
	}
}

void
FullatomDisulfideEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	core::scoring::ResSingleMinimizationData const &,
	core::scoring::ResSingleMinimizationData const &,
	core::scoring::ResPairMinimizationData & data_cache
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	// Ignore VIRTUAL residues.
	if ( rsd1.is_virtual_residue() || rsd2.is_virtual_residue() ) {
		return;
	}

	conformation::Residue const & rsdl( rsd1.seqpos() < rsd2.seqpos() ? rsd1 : rsd2 );
	conformation::Residue const & rsdu( rsd1.seqpos() < rsd2.seqpos() ? rsd2 : rsd1 );

	DisulfMinDataOP disulf_inds( new DisulfMinData( rsdl, rsdu ) );

	data_cache.set_data( core::scoring::fa_dslf_respair_data, disulf_inds );
}

void
FullatomDisulfideEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::scoring::ResSingleMinimizationData const &,
	core::scoring::ResSingleMinimizationData const &,
	core::scoring::ResPairMinimizationData const & min_data,
	pose::Pose const &,
	core::scoring::EnergyMap const & weights,
	utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
	utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	// Ignore VIRTUAL residues.
	if ( rsd1.is_virtual_residue() || rsd2.is_virtual_residue() ) {
		return;
	}

	debug_assert( utility::pointer::dynamic_pointer_cast< DisulfMinData const > ( min_data.get_data( core::scoring::fa_dslf_respair_data ) ) );

	DisulfMinData const & disulf_inds = static_cast< DisulfMinData const & > ( *min_data.get_data( core::scoring::fa_dslf_respair_data ) );

	/// this could be substantially more efficient, but there are only ever a handful of disulfides in proteins,
	/// so there's basically no point in spending time making this code faster

	//fpd old version
	if ( weights[core::scoring::dslf_ss_dst] != 0 || weights[core::scoring::dslf_cs_ang] != 0 || weights[core::scoring::dslf_ss_dih] != 0 || weights[core::scoring::dslf_ca_dih] != 0 ) {
		potential_.get_disulfide_derivatives_old(
			rsd1, rsd2, disulf_inds.res_inds( 1 ), disulf_inds.res_inds( 2 ),
			disulf_inds.res_inds(1).c_alpha_index(), weights,
			r1_atom_derivs[ disulf_inds.res_inds(1).c_alpha_index() ].f1(),
			r1_atom_derivs[ disulf_inds.res_inds(1).c_alpha_index() ].f2() );

		potential_.get_disulfide_derivatives_old(
			rsd1, rsd2, disulf_inds.res_inds( 1 ), disulf_inds.res_inds( 2 ),
			disulf_inds.res_inds(1).c_beta_index(), weights,
			r1_atom_derivs[ disulf_inds.res_inds(1).c_beta_index() ].f1(),
			r1_atom_derivs[ disulf_inds.res_inds(1).c_beta_index() ].f2() );

		potential_.get_disulfide_derivatives_old(
			rsd1, rsd2, disulf_inds.res_inds( 1 ), disulf_inds.res_inds( 2 ),
			disulf_inds.res_inds(1).disulf_atom_index(), weights,
			r1_atom_derivs[ disulf_inds.res_inds(1).disulf_atom_index() ].f1(),
			r1_atom_derivs[ disulf_inds.res_inds(1).disulf_atom_index() ].f2() );

		potential_.get_disulfide_derivatives_old(
			rsd2, rsd1, disulf_inds.res_inds( 2 ), disulf_inds.res_inds( 1 ),
			disulf_inds.res_inds(2).c_alpha_index(), weights,
			r2_atom_derivs[ disulf_inds.res_inds(2).c_alpha_index() ].f1(),
			r2_atom_derivs[ disulf_inds.res_inds(2).c_alpha_index() ].f2() );

		potential_.get_disulfide_derivatives_old(
			rsd2, rsd1, disulf_inds.res_inds( 2 ), disulf_inds.res_inds( 1 ),
			disulf_inds.res_inds(2).c_beta_index(), weights,
			r2_atom_derivs[ disulf_inds.res_inds(2).c_beta_index() ].f1(),
			r2_atom_derivs[ disulf_inds.res_inds(2).c_beta_index() ].f2() );

		potential_.get_disulfide_derivatives_old(
			rsd2, rsd1, disulf_inds.res_inds( 2 ), disulf_inds.res_inds( 1 ),
			disulf_inds.res_inds(2).disulf_atom_index(), weights,
			r2_atom_derivs[ disulf_inds.res_inds(2).disulf_atom_index() ].f1(),
			r2_atom_derivs[ disulf_inds.res_inds(2).disulf_atom_index() ].f2() );
	}

	//fpd new version
	if ( weights[ core::scoring::dslf_fa13 ] != 0 ) {
		potential_.get_disulfide_derivatives(
			rsd1, rsd2, disulf_inds.res_inds( 1 ), disulf_inds.res_inds( 2 ),
			weights, r1_atom_derivs, r2_atom_derivs );
	}
}

void
FullatomDisulfideEnergy::old_eval_atom_derivative(
	id::AtomID const & atomid,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	core::scoring::ScoreFunction const & sfxn,
	core::scoring::EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	//fpd Im not even sure if this function is used?  Anyway, it doesn't work with the new dslf_fa13 potential
	runtime_assert( !sfxn.has_nonzero_weight(core::scoring::dslf_fa13) );

	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( pose.residue( atomid.rsd() ).has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	// Ignore VIRTUAL residues.
	if ( pose.residue( atomid.rsd()).is_virtual_residue() ) {
		return;
	}

#ifdef NDEBUG
	core::scoring::disulfides::FullatomDisulfideEnergyContainerCOP dec (
		utility::pointer::static_pointer_cast< core::scoring::disulfides::FullatomDisulfideEnergyContainer const > (
			pose.energies().long_range_container( core::scoring::methods::fa_disulfide_energy )
		)
	);
#else
	core::scoring::disulfides::FullatomDisulfideEnergyContainerCOP dec (
		utility::pointer::dynamic_pointer_cast< core::scoring::disulfides::FullatomDisulfideEnergyContainer const > (
		pose.energies().long_range_container( core::scoring::methods::fa_disulfide_energy )
		)
	);
	debug_assert( dec != nullptr );
#endif

	if ( ! dec->residue_forms_disulfide( atomid.rsd() ) ) return;

	if ( dec->disulfide_atom_indices( atomid.rsd() ).atom_gets_derivatives( atomid.atomno() ) ) {
		conformation::Residue const & res( pose.residue( atomid.rsd() ));
		Vector f1( 0.0 ), f2( 0.0 );
		potential_.get_disulfide_derivatives_old(
			res,
			pose.residue( dec->other_neighbor_id( atomid.rsd()) ),
			dec->disulfide_atom_indices( atomid.rsd() ),
			dec->other_neighbor_atom_indices( atomid.rsd() ),
			atomid.atomno(),
			weights,
			f1, f2 );
		F1 += f1;
		F2 += f2;
	}
}


Real
FullatomDisulfideEnergy::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const &
) const
{
	return 0.0;
}

void
FullatomDisulfideEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}


// TwoBodyEnergy Methods:

void
FullatomDisulfideEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &sfxn,
	core::scoring::EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	// Ignore VIRTUAL residues.
	if ( rsd1.is_virtual_residue() || rsd2.is_virtual_residue() ) {
		return;
	}

	if ( ( !rsd1.type().is_disulfide_bonded() && !rsd1.type().is_sidechain_thiol() )
			|| ( !rsd2.type().is_disulfide_bonded() && !rsd2.type().is_sidechain_thiol() ) ) return;

#ifdef NDEBUG
	core::scoring::disulfides::FullatomDisulfideEnergyContainerCOP dec(
		utility::pointer::static_pointer_cast< core::scoring::disulfides::FullatomDisulfideEnergyContainer const > (
			pose.energies().long_range_container( core::scoring::methods::fa_disulfide_energy )
		)
	);
#else
	core::scoring::disulfides::FullatomDisulfideEnergyContainerCOP dec(
		utility::pointer::dynamic_pointer_cast< core::scoring::disulfides::FullatomDisulfideEnergyContainer const > (
		pose.energies().long_range_container( core::scoring::methods::fa_disulfide_energy )
		)
	);
	debug_assert( dec != nullptr );
#endif
	if ( ! dec->residue_forms_disulfide( rsd1.seqpos() ) ||
			dec->other_neighbor_id( rsd1.seqpos() ) != (Size) rsd2.seqpos() ) {
		return;
	}

	//fpd old version
	if ( sfxn.has_nonzero_weight(core::scoring::dslf_ss_dst) || sfxn.has_nonzero_weight(core::scoring::dslf_cs_ang) || sfxn.has_nonzero_weight(core::scoring::dslf_ss_dih) || sfxn.has_nonzero_weight(core::scoring::dslf_ca_dih) ) {
		Energy distance_score_this_disulfide;
		Energy csangles_score_this_disulfide;
		Energy dihedral_score_this_disulfide;
		Energy ca_dihedral_sc_this_disulf;
		bool truefalse_fa_disulf;

		potential_.score_this_disulfide_old(
			rsd1, rsd2,
			dec->disulfide_atom_indices( rsd1.seqpos() ),
			dec->other_neighbor_atom_indices( rsd1.seqpos() ), //The function change from the above line changes which index we get;
			// if we also change which rsd we use then it turns everything upside down twice and the disulfide_atom_indices inappropriately match
			distance_score_this_disulfide,
			csangles_score_this_disulfide,
			dihedral_score_this_disulfide,
			ca_dihedral_sc_this_disulf,
			truefalse_fa_disulf
		);

		// amw
		/*std::cout << "scoring dslf old: " <<std::endl;
		std::cout << emap[ core::scoring::dslf_ss_dst ] << " += " << distance_score_this_disulfide << std::endl;
		std::cout << emap[ core::scoring::dslf_cs_ang ] << " += " << csangles_score_this_disulfide << std::endl;
		std::cout << emap[ core::scoring::dslf_ss_dih ] << " += " << dihedral_score_this_disulfide << std::endl;
		std::cout << emap[ core::scoring::dslf_ca_dih ] << " += " << ca_dihedral_sc_this_disulf << std::endl;
		*/
		emap[ core::scoring::dslf_ss_dst ] += distance_score_this_disulfide;
		emap[ core::scoring::dslf_cs_ang ] += csangles_score_this_disulfide;
		emap[ core::scoring::dslf_ss_dih ] += dihedral_score_this_disulfide;
		emap[ core::scoring::dslf_ca_dih ] += ca_dihedral_sc_this_disulf;
	}

	//fpd new version
	if ( sfxn.has_nonzero_weight(core::scoring::dslf_fa13) ) {
		Energy score_i;
		potential_.score_this_disulfide(
			rsd1, rsd2,
			dec->disulfide_atom_indices( rsd1.seqpos() ), dec->other_neighbor_atom_indices( rsd1.seqpos() ),
			score_i );
		// amw

		/*std::cout << "scoring dslf new: " <<std::endl;
		std::cout << emap[ core::scoring::dslf_fa13 ] << " += " << score_i << std::endl;
		*/
		emap[ core::scoring::dslf_fa13 ] += score_i;
	}
}


bool
FullatomDisulfideEnergy::defines_intrares_energy( core::scoring::EnergyMap const & ) const
{
	return false;
}


void
FullatomDisulfideEnergy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &
) const
{
}

// LongRangeTwoBodyEnergy methods
core::scoring::methods::LongRangeEnergyType
FullatomDisulfideEnergy::long_range_type() const
{
	return core::scoring::methods::fa_disulfide_energy;
}

bool
FullatomDisulfideEnergy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size res1,
	Size res2
) const
{

	using namespace core::scoring::methods;
	if ( ! pose.energies().long_range_container( fa_disulfide_energy ) ) return false;

#ifdef NDEBUG
	core::scoring::disulfides::FullatomDisulfideEnergyContainerCOP dec(
		utility::pointer::static_pointer_cast< core::scoring::disulfides::FullatomDisulfideEnergyContainer const > (
			pose.energies().long_range_container( fa_disulfide_energy )
		)
	);
#else
	core::scoring::disulfides::FullatomDisulfideEnergyContainerCOP dec(
		utility::pointer::dynamic_pointer_cast< core::scoring::disulfides::FullatomDisulfideEnergyContainer const > (
		pose.energies().long_range_container( fa_disulfide_energy )
		)
	);
	debug_assert(dec != nullptr );
#endif
	return dec->disulfide_bonded( res1, res2 );
}

core::Size
FullatomDisulfideEnergy::version() const
{
	return 1; // Initial versioning
}

} // namespace energy_methods
} // namespace core

