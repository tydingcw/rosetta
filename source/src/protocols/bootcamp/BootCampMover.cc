// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/bootcamp/BootCampMover.cc
/// @brief The Bootcamp Mover
/// @author tydingcw (claiborne.w.tydings@vanderbilt.edu)

// Unit headers
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/bootcamp/BootCampMoverCreator.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <iostream>
#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <devel/init.hh>
#include <basic/options/option.hh>
//#include <utility/pointer/owning_ptr.hh>
//#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/random/random.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pose/variant_util.hh>
#include  <protocols/bootcamp/fold_tree_from_ss.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR( "protocols.bootcamp.BootCampMover" );

namespace protocols {
namespace bootcamp {

	/////////////////////
	/// Constructors  ///
	/////////////////////

/// @brief Default constructor
BootCampMover::BootCampMover():
	protocols::moves::Mover( BootCampMover::mover_name() )
{
    num_iterations_ = 100;
    sfxn_ = core::scoring::get_score_function();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
BootCampMover::~BootCampMover(){}

////////////////////////////////////////////////////////////////////////////////
	/// Mover Methods ///
	/////////////////////

//getter and setter for sfxn and num_iterations_
//core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
void BootCampMover::set_sfxn(core::scoring::ScoreFunctionOP sfxn) {
    sfxn_ = sfxn;
}

core::scoring::ScoreFunctionOP BootCampMover::get_sfxn() const {
    return sfxn_;
}

void BootCampMover::set_num_iterations(core::Real num_iterations) {
    num_iterations_ = num_iterations;
}
core::Size BootCampMover::get_num_iterations() const{
    return num_iterations_;
} //function is const

        /// @brief Apply the mover
void
BootCampMover::apply( core::pose::Pose& mypose){

  using namespace protocols::moves;

  //std::cout << "Hello World!" << std::endl;

  //devel::init( argc, argv );
  //utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
  //if ( filenames.size() > 0 ) {
  //    std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
  //} else {
  //    std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
  //    return 1;
  //}

  //core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
  //core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
  core::Real score = sfxn_->score( mypose );
  std::cout << "Score: " << score << std::endl;

  //modify foldtree
  auto my_tree = protocols::bootcamp::fold_tree_from_ss( mypose );
  mypose.fold_tree(my_tree);

  //chainbreak terms
  sfxn_->set_weight(core::scoring::linear_chainbreak, 1);
  core::pose::correctly_add_cutpoint_variants(mypose);

  core::Size n_resi = mypose.size();

  //core::Real gauss = numeric::random::gaussian();
  //std::cout << gauss << std::endl;

  core::Real uni = numeric::random::uniform();
  //std::cout << uni << std::endl;

  //core::Size randres = n_resi * uni + 1;
  //std::cout << randres << " randres"  << std::endl;


  MonteCarlo monte_carlo = MonteCarlo(mypose, *sfxn_, 0.6);

  //protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( mypose, true, 0 );
  //the_observer->pymol().apply( mypose);

  core::kinematics::MoveMap mm;
  mm.set_bb( true );
  mm.set_chi( true );
  core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
  core::optimization::AtomTreeMinimizer atm;

  core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( mypose );
  repack_task->restrict_to_repacking();

  core::pose::Pose copy_pose;

  bool accept;
  core::Real accept_count = 0;
  core::Real accept_ratio;
  core::Real sum_scores = 0.0;
  core::Real avg_score;

  for (core::Size i = 1; i <= num_iterations_; i++){

      core::Size randres = n_resi * uni + 1;
      core::Real pert1 = numeric::random::gaussian();
      core::Real pert2 = numeric::random::gaussian();

      core::Real orig_phi = mypose.phi( randres );
      core::Real orig_psi = mypose.psi( randres );
      mypose.set_phi( randres, orig_phi + pert1 );
      mypose.set_psi( randres, orig_psi + pert2 );

      core::pack::pack_rotamers( mypose, *sfxn_, repack_task );

      //atm.run( mypose, mm, *sfxn, min_opts );

      copy_pose = mypose;
      atm.run( copy_pose, mm, *sfxn_, min_opts );
      mypose = copy_pose;

      //Record Acceptance ratio and score
      accept = monte_carlo.boltzmann(mypose);
      if (accept) {
          accept_count++;
      }
      score = sfxn_->score( mypose );
      sum_scores += score;
      //Acceptance Ratio at Step 100: 0.81
      //Average Score at Step 100: -246.001
      if (i%100 == 0) {
          accept_ratio = accept_count/i;
          std::cout << "Acceptance Ratio at Step " << i << ": " << accept_ratio << std::endl;
          avg_score = sum_scores/i;
          std::cout << "Average Score at Step " << i << ": " << avg_score << std::endl;
      }

  }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
BootCampMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
BootCampMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap& datamap
) {

  if ( tag->hasOption("num_iterations") ) {
      num_iterations_ = tag->getOption<core::Size>("num_iterations",100);
      runtime_assert( num_iterations_ > 0 );
  }
  parse_score_function( tag, datamap );
  //parse_task_operations( tag, datamap );

}

/// @brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
void
BootCampMover::parse_score_function(
        TagCOP const tag,
        basic::datacache::DataMap const & datamap
)
{
    core::scoring::ScoreFunctionOP new_score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
    if ( new_score_function == nullptr ) return;
    set_sfxn(new_score_function);
    //sfxn_( new_score_function );
}

/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
//void
//BootCampMover::parse_task_operations(
//        TagCOP const tag,
//        basic::datacache::DataMap const & datamap
//)
//{
//    core::pack::task::TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
//    if ( new_task_factory == nullptr ) return;
//    task_factory( new_task_factory );
//}

//void BootCampMover::task_factory( core::pack::task::TaskFactoryCOP tf )
//{
//    runtime_assert( tf != nullptr );
//    task_factory_ = tf;
//}

//        void
//        BootCampMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
//            using namespace utility::tag;
//            AttributeList attlist;
//
//        }

void BootCampMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.
    attlist
    + XMLSchemaAttribute(
            "num_iterations", xsct_non_negative_integer,
            "Number of iterations for Monte Carlo Sampling" );
    rosetta_scripts::attributes_for_parse_score_function_w_description(
            attlist,
            "Use your fav score fxn");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "The Bootcamp Mover", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::fresh_instance() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::clone() const
{
	return utility::pointer::make_shared< BootCampMover >( *this );
}

std::string BootCampMover::get_name() const {
	return mover_name();
}

std::string BootCampMover::mover_name() {
	return "BootCampMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
BootCampMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

std::string
BootCampMoverCreator::keyname() const
{
	return BootCampMover::mover_name();
}

void BootCampMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
   BootCampMover::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns tydingcw as its author.
void
BootCampMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"BootCampMover", basic::citation_manager::CitedModuleType::Mover,
		"tydingcw",
		"TODO: institution",
		"claiborne.w.tydings@vanderbilt.edu",
		"Wrote the BootCampMover."
		)
	);
}


////////////////////////////////////////////////////////////////////////////////
	/// private methods ///
	///////////////////////


std::ostream &
operator<<( std::ostream & os, BootCampMover const & mover )
{
	mover.show(os);
	return os;
}


} //bootcamp
} //protocols
