// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/drug_design/bcl/BCLFragmentBaseMover.fwd.hh
/// @brief A base mover for the perturbation of BCL fragments (small molecules)
/// @author Benjamin P. Brown (benjamin.p.brown17@gmail.com)

#ifndef INCLUDED_protocols_drug_design_bcl_BCLFragmentBaseMover_fwd_hh
#define INCLUDED_protocols_drug_design_bcl_BCLFragmentBaseMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace drug_design {
namespace bcl {

class BCLFragmentBaseMover;

using BCLFragmentBaseMoverOP = utility::pointer::shared_ptr< BCLFragmentBaseMover >;
using BCLFragmentBaseMoverCOP = utility::pointer::shared_ptr< BCLFragmentBaseMover const >;

} // bcl
} //drug_design
} //protocols

#endif //INCLUDED_protocols_drug_design_bcl_BCLFragmentBaseMover_fwd_hh
