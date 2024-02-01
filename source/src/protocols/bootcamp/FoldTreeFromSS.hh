// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Utility headers
#include <core/kinematics/FoldTree.hh>
#include <protocols/loops/Loop.hh>

/// Project headers
//#include <core/types.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/moves/DsspMover.hh>

// C++ headers
#include <utility/vector1.hh>
#include <string>
#include <iostream>

namespace protocols{
namespace bootcamp{

    class FoldTreeFromSS {
    public:
        //Constructor to make ft from string
        FoldTreeFromSS( std::string const & ss_string );

        //Constructor to make ft from pose
        FoldTreeFromSS( core::pose::Pose & pose );

        core::kinematics::FoldTree const &
        fold_tree() const;

        protocols::loops::Loop const &
        loop( core::Size index ) const;

        core::Size
        loop_for_residue( core::Size seqpos ) const;

    private:
        core::kinematics::FoldTree ft_;
        utility::vector1< protocols::loops::Loop > loop_vector_;
        utility::vector1< core::Size > loop_for_residue_;
    };

/// @brief generate spans of SS, returned as vector of pairs
utility::vector1< std::pair< core::Size, core::Size > >
identify_secondary_structure_spans( std::string const & ss_string );

/// @brief generate a ft from a string, using secondary structure
core::kinematics::FoldTree fold_tree_from_dssp_string(
        std::string input_string
) ;

/// @brief generate a ft from a pose, using secondary structure
core::kinematics::FoldTree fold_tree_from_ss(core::pose::Pose & pose);

}
}
