// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/pack_rotamers.cc
/// @brief  pack rotamers module
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitue.org) to allow multi-threaded interaction graph setup.

// Unit Headers
#include <core/pack/pack_rotamers.hh>

// Package Headers
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>

#include <core/pack/annealer/AnnealerFactory.hh>
#include <core/pack/annealer/SimAnnealerBase.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>

//#include <core/kinematics/FoldTree.hh>
//#include <core/kinematics/Jump.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>


// Project Headers
#include <core/conformation/Residue.hh>
#include <utility/graph/Graph.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// util
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

// option key includes

#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <utility/thread/backwards_thread_local.hh> //For THREAD_LOCAL

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/FArray1D.hh> // AUTO IWYU For FArray1D, FArray1D<>::size_type, FArray1D::IR

using namespace ObjexxFCL;

namespace core {
namespace pack {

using core::conformation::symmetry::SymmetryInfoCOP;

static basic::Tracer tt( "core.pack.pack_rotamers", basic::t_info );

/// @brief The entry point to calling the packer, which optimizes side-chain identities and conformations
/// (if you're designing) or conformations alone (if you're predicting structures).
/// @details Wraps the two very distinct and separate stages of rotamer packing,
/// which are factored so that they may be called asynchronously.
/// Use this wrapper as a base model for higher-level packing routines (such as pack_rotamers_loop).
/// @note In multi-threaded builds, this function takes an extra parameter for
/// the number of threads to request, for parallel interaction graph precomputation.
void
pack_rotamers(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task
) {
	using namespace interaction_graph;
	using namespace rotamer_set;

	task::PackerTaskCOP working_task = nullptr;
	rotamer_set::RotamerSetsOP rotsets = nullptr;

	//fpd safety check for symmetry
	//replace this with RotSetsFactory
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		working_task = make_new_symmetric_PackerTask_by_requested_method(pose, task);
		rotsets = utility::pointer::make_shared< rotamer_set::symmetry::SymmetricRotamerSets >();
	} else {
		working_task = task;
		rotsets = utility::pointer::make_shared< rotamer_set::RotamerSets >();
	}

	PROF_START( basic::PACK_ROTAMERS );

	pack_scorefxn_pose_handshake( pose, scfxn);

	AnnealableGraphBaseOP ig = nullptr;

	pack_rotamers_setup( pose, scfxn, task, rotsets, ig );
	pack_rotamers_run( pose, task, rotsets, ig );
	pack_rotamers_cleanup( pose, ig );

	// rescore here to make the state of the Energies good.
	scfxn( pose );

	PROF_STOP ( basic::PACK_ROTAMERS );
}

/// @details run the FixbbSimAnnealer multiple times using the same InteractionGraph, storing the results.
/// @note In multi-threaded builds, this function takes an extra parameter for
/// the number of threads to request, for parallel interaction graph precomputation.
void
pack_rotamers_loop(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	Size const nloop
) {
	utility::vector1< std::pair< Real, std::string > > results;
	pack_rotamers_loop( pose, scfxn, task, nloop, results );
}

/// @details run the FixbbSimAnnealer multiple times using the same InteractionGraph, storing the results
/// @note In multi-threaded builds, this function takes an extra parameter for
/// the number of threads to request, for parallel interaction graph precomputation.
void
pack_rotamers_loop(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	Size const nloop,
	utility::vector1< std::pair< Real, std::string > > & results
) {
	utility::vector1< pose::PoseOP > pose_list;
	pack_rotamers_loop( pose, scfxn, task, nloop, results, pose_list );
}

/// @brief Run the FixbbSimAnnealer multiple times using the same InteractionGraph, storing the results.
/// @note In multi-threaded builds, this function takes an extra parameter for
/// the number of threads to request, for parallel interaction graph precomputation.
void
pack_rotamers_loop(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	Size const nloop,
	utility::vector1< std::pair< Real, std::string > > & results,
	utility::vector1< pose::PoseOP > & pose_list
) {
	using namespace ObjexxFCL::format;
	using namespace interaction_graph;
	using namespace rotamer_set;

	rotamer_set::RotamerSetsOP rotsets( new rotamer_set::RotamerSets() );
	AnnealableGraphBaseOP ig = nullptr;
	pack_rotamers_setup( pose, scfxn, task, rotsets, ig );

	Real best_bestenergy( 0.0 );
	pose::Pose best_pose;
	best_pose = pose;

	for ( Size run(1); run <= nloop; ++run ) {

		Real bestenergy( pack_rotamers_run( pose, task, rotsets, ig ) );

		Real const final_score( scfxn( pose ) );

		// show the resulting sequence
		std::string final_seq;
		for ( Size i=1; i<= pose.size(); ++i ) {
			if ( task->design_residue(i) ) final_seq+= pose.residue(i).name1();
		}
		if ( final_seq.size() == 0 ) final_seq = "-";

		tt << "packloop: " << I(4,run) << " rescore-deltaE: " << F(9,3,final_score - bestenergy ) <<
			" seq: " << final_seq << " simannealerE: " << bestenergy << " rescoreE: " << final_score << std::endl;

		results.push_back( std::make_pair( final_score, pose.sequence() ) );
		pose_list.push_back( utility::pointer::make_shared< pose::Pose >( pose ) ); //saves a copy.
		if ( run == 1 || bestenergy < best_bestenergy ) {
			best_pose = pose;
			best_bestenergy = bestenergy;
		}
	}
	pose = best_pose;
	pack_rotamers_cleanup( pose, ig );
}

/// @brief Get rotamers, compute energies.
/// @note In multi-threaded builds, this function takes an extra parameter for
/// the number of threads to request, for parallel interaction graph precomputation.
void
pack_rotamers_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP rotsets,
	interaction_graph::AnnealableGraphBaseOP & ig,
	core::Size const nloop //how many packing runs will use this IG? This totally optional parameter helps the factory decide which IG to use if nothing else is specified
) {
	using namespace interaction_graph;

	pack_scorefxn_pose_handshake( pose, scfxn );

	pose.update_residue_neighbors();

	scfxn.setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );

	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scfxn, task );

	rotsets->set_task( task );
	rotsets->initialize_pose_for_rotsets_creation(pose);
	rotsets->build_rotamers( pose, scfxn, packer_neighbor_graph );
	scfxn.setup_for_packing_with_rotsets( pose, rotsets );
	rotsets->prepare_sets_for_packing( pose, scfxn );

	if ( basic::options::option[ basic::options::OptionKeys::packing::dump_rotamer_sets ] ) { // hacking
#ifndef MULTI_THREADED
		static THREAD_LOCAL int counter(0);
		++counter;
		std::string const filename( "rotset"+lead_zero_string_of( counter,4 )+".pdb" );
		tt << "dump rotsets: " << filename << std::endl;
		rotsets->dump_pdb( pose, filename );
#else
		utility_exit_with_message("Error: the -packing:dump_rotamer_sets flag cannot be used in the multi-threaded build.");
#endif
	}

	tt << "built " << rotsets->nrotamers() << " rotamers at "
		<< rotsets->nmoltenres() << " positions." << std::endl;

#ifdef MULTI_THREADED
	core::Size const threads_to_request( task->ig_threads_to_request() );
	tt << "Requesting ";
	if ( threads_to_request == 0 ) {
		tt << "all available";
	} else {
		tt << threads_to_request;
	}
	tt << " threads for interaction graph computation." << std::endl;
#endif
	ig = InteractionGraphFactory::create_and_initialize_annealing_graph(*task, *rotsets, pose, scfxn, packer_neighbor_graph, nloop );
}

/// @brief PyRosetta compatible version.
interaction_graph::AnnealableGraphBaseOP
pack_rotamers_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP rotsets
) {
	interaction_graph::AnnealableGraphBaseOP ig;
	pack_rotamers_setup(pose, scfxn, task, rotsets, ig);
	return ig;
}

// @begin pack_rotamers_run
// @details as simple as possible -- runs simulated annealing and then places
// the optimal rotamers onto the backbone of the input pose.
Real
pack_rotamers_run(
	pose::Pose & pose,
	task::PackerTaskCOP task,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::AnnealableGraphBaseOP ig,
	utility::vector0< int > rot_to_pack, // defaults to an empty vector (no effect)
	annealer::AnnealerObserverOP observer //defaults to nullptr
)
{
	using namespace ObjexxFCL::format;

	FArray1D_int bestrotamer_at_seqpos( pose.size() );
	core::PackerEnergy bestenergy( 0.0 );

	pack_rotamers_run( pose, task, rotsets, ig, rot_to_pack, observer, bestrotamer_at_seqpos, bestenergy );

	// place new rotamers on input pose
	for ( uint ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
		uint iiresid = rotsets->moltenres_2_resid( ii );
		uint iibestrot = rotsets->rotid_on_moltenresidue( bestrotamer_at_seqpos( iiresid ) );
		conformation::ResidueCOP bestrot( rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( iibestrot ) );

		conformation::ResidueOP newresidue( bestrot->create_residue() );
		pose.replace_residue ( iiresid, *newresidue, false );
	}

	return bestenergy;
}

/// @brief Runs simulated annealing and returns the best identified rotamers and energies.
void
pack_rotamers_run(
	pose::Pose const & pose,
	task::PackerTaskCOP task,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::AnnealableGraphBaseOP ig,
	utility::vector0< int > rot_to_pack,
	annealer::AnnealerObserverOP observer,
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy
)
{
	using namespace annealer;

	bool start_with_current = false;
	FArray1D_int current_rot_index( pose.size(), 0 );
	bool calc_rot_freq = false;
	FArray1D< core::PackerEnergy > rot_freq( ig->get_num_total_states(), 0.0 );

	// too many parameters for annealer's constructor! should replace with a task only;
	// the annealer should then provide read access to the data its collected after
	// annealing has completed. some data it won't bother collecting because the task
	// did not instruct it to.

	/// Parameters passed by reference in task's constructor to which it writes at the
	/// completion of sim annealing.

	SimAnnealerBaseOP annealer = AnnealerFactory::create_annealer(
		pose, task, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current,
		ig, rotsets, current_rot_index, calc_rot_freq, rot_freq );

	// Following additional initialization would be cleaner if above suggestion
	// regarding the removal of "too many parameters" is implemented properly.
	if ( task->low_temp()  > 0.0 ) annealer->set_lowtemp(  task->low_temp()  );
	if ( task->high_temp() > 0.0 ) annealer->set_hightemp( task->high_temp() );
	annealer->set_disallow_quench( task->disallow_quench() );
	annealer->set_annealer_observer( observer );

	PROF_START( basic::SIMANNEALING );
	annealer->run();
	PROF_STOP( basic::SIMANNEALING );
}

/// @brief Provide the opportunity to clean up cached data from the pose or scorefunction after packing.
/// @details This should be called after pack_rotamers_run.  It is called from the pack_rotamers() function.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
pack_rotamers_cleanup(
	core::pose::Pose & pose,
	core::pack::interaction_graph::AnnealableGraphBaseOP annealable_graph
) {
	annealable_graph->clean_up_after_packing(pose);
}

} // namespace pack
} // namespace core
