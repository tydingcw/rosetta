// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SymmetryClaim
/// @brief Claim for Symmetry (includes SymmData)
/// @author Justin Porter, Tatjana Braun


#ifndef INCLUDED_protocols_topology_broker_claims_SymmetryClaim_hh
#define INCLUDED_protocols_topology_broker_claims_SymmetryClaim_hh


// Unit Headers
#include <protocols/topology_broker/claims/SymmetryClaim.fwd.hh>
#include <protocols/topology_broker/claims/DofClaim.hh>


// Package Headers
#include <protocols/topology_broker/TopologyClaimer.fwd.hh> //for printing

// Project Headers
#include <core/conformation/symmetry/SymmData.fwd.hh>

// ObjexxFCL Headers

// Utility headers

//// C++ headers
#include <string>
#include <iosfwd>

#include <ostream> // AUTO IWYU For operator<<, basic_ostream

namespace protocols {
namespace topology_broker {
namespace claims {

//this is a bit different then the other claims.
//overlapping sequences don't make sense I think.
// this is merely used to arrange multiple patches so that everybody has a residue number.
// pos defines the starting residue of the patch.
// length the length.
// the initial claim goes out with pos=0 if you allow this sequence to move.
// length specifies the number of residues.
// pos=0 claims will be assigned a position from the broker.
class SymmetryClaim : public DofClaim {
public:
	SymmetryClaim( TopologyClaimerAP tc, core::conformation::symmetry::SymmDataOP symmdata,
		std::string const& label, ClaimRight right )  :
		DofClaim( tc, right ),
		label_( label ),
		symm_data_( symmdata )
	{}

	DofClaimOP clone() const override {
		return utility::pointer::make_shared< SymmetryClaim >( *this );
	}

	std::string const& label() const {
		return label_;
	}

	core::conformation::symmetry::SymmDataOP get_symm_data(){
		return symm_data_;
	}

	std::string str_type() const override {
		return "SYMMETRY";
	}

	void show( std::ostream& os ) const override {
		os << " with label: " << label();
	}

	//    virtual std::string to_string() const {
	//        std::ostringstream str_stream;
	//        str_stream << "(SymmetryClaim; owner, " << owner()->type() //<< "; label, " << label()
	//            << "; symmetry name, " << symm_data_->get_symmetry_name() << ")" ;
	//        return str_stream.str();
	//    }

private:
	std::string label_;
	core::conformation::symmetry::SymmDataOP symm_data_;
};

} // claims
} // topology_broker
} // protocols
#endif
