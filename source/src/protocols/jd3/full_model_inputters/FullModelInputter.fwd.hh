// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/full_model_inputters/FullModelInputter.fwd.hh
/// @brief  Forward declaration of the %FullModelInputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_full_model_inputters_FullModelInputter_FWD_HH
#define INCLUDED_protocols_jd3_full_model_inputters_FullModelInputter_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd3 {
namespace full_model_inputters {

class FullModelInputter;

typedef utility::pointer::shared_ptr< FullModelInputter > FullModelInputterOP;
typedef utility::pointer::shared_ptr< FullModelInputter const > FullModelInputterCOP;

} // namespace full_model_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_FullModelInputter_FWD_HH
