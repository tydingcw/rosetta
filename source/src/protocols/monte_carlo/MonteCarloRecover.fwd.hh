// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_MonteCarloRecover_fwd_hh
#define INCLUDED_protocols_simple_moves_MonteCarloRecover_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace monte_carlo {

class MonteCarloRecover;
typedef utility::pointer::shared_ptr< MonteCarloRecover >  MonteCarloRecoverOP;
typedef utility::pointer::shared_ptr< MonteCarloRecover const >  MonteCarloRecoverCOP;


} // namespace monte_carlo
} // namespace protocols

#endif
