// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/Residue.functions.hh
/// @brief  Header file for non-member functions that operate on Residue objects
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_conformation_Residue_functions_hh
#define INCLUDED_core_conformation_Residue_functions_hh

// Package headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/PartialAtomID.fwd.hh>

// C++ headers
#include <set>


namespace core {
namespace conformation {

void
idealize_hydrogens(
	Residue & res,
	Conformation const & conf
);

/// @brief rotamer chi-update from coords useful for building rotamers from coordinates
void set_chi_according_to_coordinates(
	conformation::Residue & rotamer
);

/// @brief For a given mainchain torsion, report the partial atom IDs that
/// define it, inserting these partial atom IDs into a std::set (to avoid redundancy
/// perhaps with other torsions that also depend on the same atoms).
void
insert_partial_atom_ids_for_mainchain_torsion(
	Residue const & res,
	Size const mainchain_torsion,
	std::set< id::PartialAtomID > & atoms
);

}
}

#endif
