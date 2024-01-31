// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/bootcamp/BootCampMover.fwd.hh
/// @brief The Bootcamp Mover
/// @author tydingcw (claiborne.w.tydings@vanderbilt.edu)

#ifndef INCLUDED_protocols_bootcamp_BootCampMover_fwd_hh
#define INCLUDED_protocols_bootcamp_BootCampMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace bootcamp {

class BootCampMover;

using BootCampMoverOP = utility::pointer::shared_ptr< BootCampMover >;
using BootCampMoverCOP = utility::pointer::shared_ptr< BootCampMover const >;

} //bootcamp
} //protocols

#endif //INCLUDED_protocols_bootcamp_BootCampMover_fwd_hh
