// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/GaussianOverlapEnergy.hh
/// @brief  energy of overlaps
/// @author  ben


#ifndef INCLUDED_core_energy_methods_GaussianOverlapEnergy_hh
#define INCLUDED_core_energy_methods_GaussianOverlapEnergy_hh

// Unit Headers
#include <core/energy_methods/GaussianOverlapEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <numeric/interpolation/spline/Interpolator.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>

// Utility headers


namespace core {
namespace energy_methods {



class GaussianOverlapEnergy : public core::scoring::methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentTwoBodyEnergy  parent;
public:


	GaussianOverlapEnergy();
	~GaussianOverlapEnergy() override;


	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;


	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & scorefxn,
		core::scoring::EnergyMap & emap
	) const override;


	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;


	Distance
	atomic_interaction_cutoff() const override;

	/// @details non-virtual accessor for speed
	Distance
	interaction_cutoff() const;


	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const override;

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override {
		return true;
	}

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap & emap
	) const override;


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	numeric::interpolation::spline::InterpolatorOP interp_;
	core::Size version() const override;

};


}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
