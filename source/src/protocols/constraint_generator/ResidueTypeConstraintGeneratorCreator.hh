// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/ResidueTypeConstraintGeneratorCreator.hh
/// @brief Generates residue type constraints for a set of residues from the current pose
/// @author Sharon Guffy (guffy@email.unc.edu)

#ifndef INCLUDED_protocols_constraint_generator_ResidueTypeConstraintGeneratorCreator_hh
#define INCLUDED_protocols_constraint_generator_ResidueTypeConstraintGeneratorCreator_hh

// Unit headers
#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace constraint_generator {

class ResidueTypeConstraintGeneratorCreator : public protocols::constraint_generator::ConstraintGeneratorCreator {
public:
	protocols::constraint_generator::ConstraintGeneratorOP
	create_constraint_generator() const override;

	std::string
	keyname() const override;

	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols/constraint_generator_ResidueTypeConstraintGenerator_fwd_hh
