// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/energy_methods/PointWaterEnergy.hh
/// @brief  Statistical point water energy function
/// @author Frank DiMaio


#ifndef INCLUDED_core_energy_methods_PointWaterEnergy_fwd_hh
#define INCLUDED_core_energy_methods_PointWaterEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace energy_methods {


class PointWaterEnergy;

typedef utility::pointer::shared_ptr< PointWaterEnergy > PointWaterEnergyOP;
typedef utility::pointer::shared_ptr< PointWaterEnergy const > PointWaterEnergyCOP;

}
}

#endif
