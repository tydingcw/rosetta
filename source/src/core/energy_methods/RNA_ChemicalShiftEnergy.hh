// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/RNA_ChemicalShiftEnergy.hh
/// @brief  Energy score based on the agreement between experimentally determined and theoretically calculated NMR chemical_shift
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)


#ifndef INCLUDED_core_scoring_rna_chemical_shift_RNA_ChemicalShiftEnergy_HH
#define INCLUDED_core_scoring_rna_chemical_shift_RNA_ChemicalShiftEnergy_HH

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>


// Utility headers


namespace core {
namespace energy_methods {


/////////////////////////////////////////////////////////////////////////////
class RNA_ChemicalShiftEnergy : public core::scoring::methods::WholeStructureEnergy{

public:
	typedef core::scoring::methods::WholeStructureEnergy  parent;

public:

	RNA_ChemicalShiftEnergy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	setup_for_scoring( pose::Pose &, core::scoring::ScoreFunction const & ) const override;

	void
	setup_for_derivatives( pose::Pose &, core::scoring::ScoreFunction const & ) const override;

	/////////////////////////////////
	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;

	/////////////////////////////////
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2 ) const override;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const override {}

private:

	core::scoring::rna::chemical_shift::RNA_ChemicalShiftPotential rna_chemical_shift_potential_; //Can I make this not mutable?

	core::Size version() const override;

};


} //scoring
} //core

#endif // INCLUDED_core_energy_methods_RG_Energy_RNA_HH
