// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/buns/BuriedUnsatHbondFilter2.hh
/// @brief
/// @author Kevin Houlihan (khouli@unc.edu)

#ifndef INCLUDED_protocols_buns_BuriedUnsatHbondFilter2_hh
#define INCLUDED_protocols_buns_BuriedUnsatHbondFilter2_hh

#include <protocols/buns/BuriedUnsatHbondFilter2.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/buns/BuriedUnsatisfiedPolarsCalculator2.fwd.hh>

namespace protocols {
namespace buns {


/// @brief filters based on an upper bound # of buried unsatisfied polar residues
class BuriedUnsatHbondFilter2 : public protocols::filters::Filter
{

public:
	// @brief default constructor
	BuriedUnsatHbondFilter2();
	// @brief constructor with options
	BuriedUnsatHbondFilter2( core::Size const upper_threshold, core::Size const jump_num );
	// @brief copy constructor
	BuriedUnsatHbondFilter2( BuriedUnsatHbondFilter2 const & rval );
	// @brief make clone
	protocols::filters::FilterOP clone() const override { return utility::pointer::make_shared< BuriedUnsatHbondFilter2 >( *this ); }
	// @brief make fresh instance
	protocols::filters::FilterOP fresh_instance() const override { return utility::pointer::make_shared< BuriedUnsatHbondFilter2 >(); }
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Size compute( core::pose::Pose const & pose ) const;

	~BuriedUnsatHbondFilter2() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;
	void task_factory( core::pack::task::TaskFactoryOP tf );
	core::pack::task::TaskFactoryOP task_factory() const;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::scoring::ScoreFunctionCOP sfxn_;
	core::Size upper_threshold_;
	core::Size jump_num_;
	core::pack::task::TaskFactoryOP task_factory_; // dflt NULL; only residues defined as packable by the taskoperations will be tested for burial
	//protocols::simple_pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator2OP calc_;
	protocols::buns::BuriedUnsatisfiedPolarsCalculator2OP calc_;
	// std::string calc_name_;
};

} //buns
} //protocols
#endif
