// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/count_pair/CountPairAll.hh
/// @brief  Count pair for residues where all atom pairs should be counted.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_etable_count_pair_CountPairAll_hh
#define INCLUDED_core_scoring_etable_count_pair_CountPairAll_hh

// Unit Headers
#include <core/scoring/etable/count_pair/CountPairAll.fwd.hh>

// Package Headers
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

class CountPairAll : public CountPairFunction
{
public:
public:
	typedef CountPairFunction parent;

public:
	CountPairAll() {} // inlined when declared on the stack
	~CountPairAll() override {} // inlined when declared on the stack

	/// @brief function required by templated functions in atom_pair_energy_inline
	inline
	bool
	operator () (
		int const /*at1*/,
		int const /*at2*/,
		Real & /*weight*/,
		Size & /*path_dist*/
	) const
	{
		return true ; //compiler, please remove an if(true) test
	}

	bool
	count(
		int const at1,
		int const at2,
		Real & w,
		Size & path_dist
	) const override;

	/// Type resolution functions
	void
	residue_atom_pair_energy(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const override;

	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const override;

	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const override;

	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const override;


	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const override;


	/// Type resolution functions
	void
	residue_atom_pair_energy(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const override;

	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const override;

	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const override;

	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const override;


	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const override;

};

} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

#endif
