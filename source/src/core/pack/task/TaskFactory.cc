// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/TaskFactory.hh
/// @brief  Task class to describe packer's behavior header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- added support for PackerPalettes, multithreading.

#include <utility/exit.hh>

// Unit Headers
#include <core/pack/task/TaskFactory.hh>

// Package Headers
#include <core/pack/task/PackerTask_.hh> // only place in all of mini where this gets #included (except PackerTask_.cc)
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/palette/PackerPalette.hh>
#include <core/pack/palette/PackerPaletteFactory.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/multithreading.OptionKeys.gen.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/list.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.pack.task.TaskFactory" );

namespace core {
namespace pack {
namespace task {

/// @brief Default constructor.
/// @details Defaults to storing a DefaultPackerPalette; the operations_ list starts out empty.
TaskFactory::TaskFactory() :
	parent()
{}

/// @brief Copy constructor.
///
TaskFactory::TaskFactory( TaskFactory const & src)
:
	parent()
{
	if ( src.packer_palette_ != nullptr ) {
		packer_palette_ = src.packer_palette_->clone();
	} else {
		packer_palette_ = nullptr;
	}
	copy_operations( src );
}

TaskFactoryOP TaskFactory::clone() const
{
	return utility::pointer::make_shared< TaskFactory >(*this);
}

TaskFactory::~TaskFactory() = default;

TaskFactory &
TaskFactory::operator = ( TaskFactory const & rhs )
{
	copy_operations( rhs );
	return *this;
}


/// @brief Apply each of the TaskOperations in the TaskFactory list to the PackerTask to set it up.
/// @details Must be called AFTER PackerTask initialization.  An uninitialized PackerTask
/// that is modified will throw an error.
void
TaskFactory::modify_task( core::pose::Pose const & pose, PackerTaskOP task ) const
{
	runtime_assert_string_msg( task != nullptr, "Error in core::pack::task::TaskFactory::modify_task(): a null pointer to the PackerTask was provided." );
	runtime_assert_string_msg( task->is_initialized(), "Error in core::pack::task::TaskFactory::modify_task(): a PackerTask was provided that has not yet been initialized." );
	for ( TaskOperationOP const & taskop : *this ) {
		taskop->apply( pose, *task );
	}
}

// Non static version.
PackerTaskOP
TaskFactory::create_task_and_apply_taskoperations( pose::Pose const & pose, core::Size const ig_threads_to_request /*=0*/ ) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// If a PackerPalette has not been set, use the default one.
	core::pack::palette::PackerPaletteOP packer_palette( packer_palette_ );
	if ( packer_palette == nullptr ) {
		TR.Debug << "No PackerPalette has been set up for this TaskFactory.  Getting the global default PackerPalette." << std::endl;
		packer_palette = core::pack::palette::PackerPaletteFactory::get_instance()->create_packer_palette_from_global_defaults();
	}

	PackerTaskOP task( utility::pointer::make_shared< PackerTask_ >( pose, packer_palette, (ig_threads_to_request == 0 ? option[multithreading::interaction_graph_threads]() : ig_threads_to_request ) ) );
	modify_task( pose, task );
	if ( TR.Trace.visible() ) {
		TR.Trace << "PackerTask created by TaskFactory:" << std::endl;
		TR.Trace << *task << std::endl;
		TR.Trace << "Attached a " << packer_palette->name() << " to the PackerTask." << std::endl;
	}
	return task;
}

// clones the input task, and pushes it back into the list
void
TaskFactory::push_back( TaskOperationCOP taskop )
{
	operations_.push_back( taskop->clone() );
}

/// @brief Clones the input PackerPalette, and sets it as the PackerPalette for this TaskFactory.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
TaskFactory::set_packer_palette(
	PackerPaletteCOP packer_palette_in
) {
	runtime_assert_string_msg( packer_palette_in, "Error in core::pack::task::TaskFactory::set_packer_palette(): a null pointer to the PackerPalette was provided." );
	packer_palette_ = packer_palette_in->clone();
}

TaskFactory::const_iterator
TaskFactory::begin() const
{
	return operations_.begin();
}

TaskFactory::const_iterator
TaskFactory::end() const
{
	return operations_.end();
}

/// @brief Empties the list of TaskOperations and clears the PackerPalette, replacing it with a new DefaultPackerPalette.
///
void
TaskFactory::clear()
{
	packer_palette_ = nullptr;
	operations_.clear();
}

/// @brief How many TaskOperations does this TaskFactory contain?
core::Size
TaskFactory::n_task_operations() const {
	return operations_.size();
}

/// @brief Does this TaskFactory contain TaskOperations?
/// @details Returns true if and only if the number of task operations is
/// greater than or equal to one.
bool
TaskFactory::has_task_operations() const {
	return !(operations_.empty());
}

/// @brief Static construction of a task
/// @details Returns a new PackerTask with NO TaskOperations, and a default PackerPalette applied.
PackerTaskOP
TaskFactory::create_packer_task(
	pose::Pose const & pose,
	core::Size const ig_threads_to_request /*=0*/
)
{
	return utility::pointer::make_shared< PackerTask_ >( pose, ( ig_threads_to_request == 0 ? basic::options::option[basic::options::OptionKeys::multithreading::interaction_graph_threads]() : ig_threads_to_request ) );
}

/// @brief Static construction of a task
/// @details Returns a new PackerTask with NO TaskOperations, and with a supplied
/// PackerPalette (used directly -- not cloned on input).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
PackerTaskOP
TaskFactory::create_packer_task(
	pose::Pose const & pose,
	core::pack::palette::PackerPaletteCOP palette,
	core::Size const ig_threads_to_request /*=0*/
) {
	return utility::pointer::make_shared< PackerTask_ >( pose, palette, ( ig_threads_to_request == 0 ? basic::options::option[basic::options::OptionKeys::multithreading::interaction_graph_threads]() : ig_threads_to_request ) );
}

void
TaskFactory::copy_operations( TaskFactory const & src )
{
	for ( auto const & taskop_iter : src ) {
		operations_.push_back( taskop_iter->clone() );
	}
}


core::Size
TaskFactory::size() const
{
	core::Size const size( operations_.size() );
	return( size );
}

} //namespace task
} //namespace pack
} //namespace core

#ifdef SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::task::TaskFactory::save( Archive & arc ) const {
	arc( CEREAL_NVP( packer_palette_) );
	arc( CEREAL_NVP( operations_ ) ); //OperationList
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::task::TaskFactory::load( Archive & arc ) {
	arc( packer_palette_ );
	arc( operations_ ); //OperationList
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::task::TaskFactory );
CEREAL_REGISTER_TYPE( core::pack::task::TaskFactory )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_TaskFactory )
#endif // SERIALIZATION
