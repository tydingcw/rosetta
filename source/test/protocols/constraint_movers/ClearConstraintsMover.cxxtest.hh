// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_movers/ClearConstraintsMover.cxxtest.hh
/// @brief  test for ClearConstraintsMover mover
/// @author Andrew Watkins

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers

#include <protocols/constraint_movers/ClearConstraintsMover.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstantConstraint.hh>
#include <core/scoring/func/ConstantFunc.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh> // AUTO IWYU For Pose, Pose::ConstraintSetCOP

static basic::Tracer TR("protocols.constraint_movers.ClearConstraintsMover.cxxtest.hh");

// --------------- Test Class --------------- //

class ClearConstraintsMoverTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_it_clears_constraints() {
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;

		core::pose::Pose pose;
		pose.add_constraint( utility::pointer::make_shared< ConstantConstraint >( utility::pointer::make_shared< ConstantFunc >( 5 ) ) );
		// Pose has one constraint

		protocols::constraint_movers::ClearConstraintsMover ccm;
		ccm.apply( pose );
		TS_ASSERT( pose.constraint_set()->is_empty() );
	}

};
