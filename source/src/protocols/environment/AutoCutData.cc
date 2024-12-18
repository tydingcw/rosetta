// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/environment/AutoCutData.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/AutoCutData.hh>

// Package headers

// Project headers

// Utility Headers
#include <utility/excn/Exceptions.hh>

// Basic Headers
#include <basic/Tracer.hh>

#include <protocols/environment/AutoCutDataCreator.hh> // AUTO IWYU For AutoCutDataCreator

// ObjexxFCL Headers

static basic::Tracer tr( "protocols.environment.AutoCutData", basic::t_info );

static std::string TYPE_NAME() { return "AutoCutData"; }

namespace protocols {
namespace environment {

basic::datacache::WriteableCacheableDataOP
AutoCutDataCreator::create_data( std::istream &in ) const {
	return utility::pointer::make_shared< AutoCutData >( in );
}

std::string AutoCutDataCreator::keyname() const{
	return TYPE_NAME();
}

AutoCutData::AutoCutData( std::istream &in ) :
	Parent()
{
	std::string token;

	in >> token;
	if ( token != "HASH" ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "AutoCutData tried to read an improperly formatted SilentFile remark." );
	}

	in >> hash_;

	in >> token;
	if ( token != "CUTS" ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "AutoCutData tried to read an improperly formatted SilentFile remark." );
	}

	while ( in.good() ) {
		core::Size cut;
		in >> cut;
		auto_cuts_.push_back( cut );
	}

	std::sort( auto_cuts_.begin(), auto_cuts_.end() );
}


AutoCutData::AutoCutData( core::Size const& hash,
	std::set< core::Size > const& cuts ) :
	Parent(),
	auto_cuts_(),
	hash_( hash )
{
	std::copy( cuts.begin(), cuts.end(), std::back_inserter( auto_cuts_ ) );
}

void AutoCutData::write( std::ostream &out ) const {

	out << TYPE_NAME() ;
	out << " HASH ";
	out << hash() << " ";

	out << "CUTS ";
	for ( core::Size const cut : cuts() ) {
		out << cut << " ";
	}
}

basic::datacache::CacheableDataOP
AutoCutData::clone() const {
	return utility::pointer::make_shared< AutoCutData >( *this );
}

std::string AutoCutData::datatype() const {
	return TYPE_NAME();
}


} // environment
} // protocols
