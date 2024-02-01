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


#include <protocols/bootcamp/FoldTreeFromSS.hh>

// Utility headers
#include <core/kinematics/FoldTree.hh>

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

//Constructor to make ft from string
FoldTreeFromSS::FoldTreeFromSS( std::string const & ss_string ){
        fold_tree_from_dssp_string(ss_string); //returns ft
}

//Constructor to make ft from pose
FoldTreeFromSS::FoldTreeFromSS( core::pose::Pose & pose ){
        fold_tree_from_ss( pose ); // returns ft
}

core::kinematics::FoldTree const &
FoldTreeFromSS::fold_tree() const{
    return ft_;
}

protocols::loops::Loop const &
FoldTreeFromSS::loop( core::Size index ) const{
    return loop_vector_[index];
}

core::Size
FoldTreeFromSS::loop_for_residue( core::Size seqpos ) const{
    return loop_for_residue_[seqpos];
}

utility::vector1< std::pair< core::Size, core::Size > >
identify_secondary_structure_spans( std::string const & ss_string )
{
    utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
    core::Size strand_start = -1;
    for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
        if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
            if ( int( strand_start ) == -1 ) {
                strand_start = ii;
            } else if ( ss_string[ii] != ss_string[strand_start] ) {
                ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
                strand_start = ii;
            }
        } else {
            if ( int( strand_start ) != -1 ) {
                ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
                strand_start = -1;
            }
        }
    }
    if ( int( strand_start ) != -1 ) {
        // last residue was part of a ss-eleemnt
        ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
    }
    for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
        std::cout << "SS Element " << ii << " from residue "
                  << ss_boundaries[ ii ].first << " to "
                  << ss_boundaries[ ii ].second << std::endl;
    }
    return ss_boundaries;
}

//core::kinematics::FoldTree
void FoldTreeFromSS::fold_tree_from_dssp_string(
        const std::string& input_string
        //utility::vector1< std::pair< core::Size, core::Size > > &input_vect
) {
    auto input_vect = identify_secondary_structure_spans(input_string);
    std::cout << "New string length " << input_string.size() << std::endl;
    std::cout << "New vector length " << input_vect.size() << std::endl;

    core::kinematics::FoldTree ft;
    //core::Size edges = 4*(input_vect.size()-2);
    //core::Size jumps = 2*(input_vect.size()-2);

    //making ss_stretch, the segments that will form fold tree regions
    utility::vector1< std::pair< core::Size, core::Size > > ss_stretch;
    core::Size prior_end = input_vect[1].second;
    ss_stretch.push_back( std::make_pair( 1, prior_end));
    for (size_t i = 2; i <= input_vect.size() - 1; i++){
        core::Size curr_start = input_vect[i].first;
        core::Size curr_end = input_vect[i].second;
        if (curr_end - curr_start >=3){ //if not true, deal with this later
            if (curr_start - prior_end >= 5){ //Looking at size 3 (start-1, end+1)
                ss_stretch.push_back( std::make_pair( prior_end+1, curr_start-1) );
                ss_stretch.push_back( std::make_pair( curr_start, curr_end) );
                prior_end = curr_end;
            } else {
               //prior loop not of adequete size, merge together
                ss_stretch.push_back( std::make_pair( prior_end+1, curr_end) );
                prior_end = curr_end;
            }
        } else {
            std::cout << "Warning: did not expect to be here" << std::endl;
        }
    }
    //Now need to deal with the last ss element pair (input_string.size())
    core::Size final_start = input_vect[input_vect.size()].first;
    core::Size final_end = input_string.size();
    if (final_end - final_start >=3){ //if not true, deal with this later
        if (final_start - prior_end >= 5){ //Looking at size 3 (start-1, end+1)
            ss_stretch.push_back( std::make_pair( prior_end+1, final_start-1) );
            ss_stretch.push_back( std::make_pair( final_start, final_end) );
            //prior_end = curr_end;
        } else {
            //prior loop not of adequete size
            ss_stretch.push_back( std::make_pair( prior_end+1, final_end) );
            //prior_end = curr_end;
        }
    } else {
        std::cout << "Error: did not expect to be here (final)" << std::endl;
    }

    //Lets print to thw screen to double check output
    std::cout << "Size: " << final_end << std::endl;
    for (size_t i = 1; i <= ss_stretch.size(); i++){
        std::cout << "pair " << ss_stretch[i].first << " " << ss_stretch[i].second << std::endl;
    }

    //Making the foldtree start
    std::pair< core::Size, core::Size > node1 = ss_stretch[1];
    core::Size center = (node1.first + node1.second)/2;
    std::cout << "New edge " << center << " " << 1 << std::endl;
    std::cout << "New edge " << center << " " << node1.second << std::endl;
    ft.add_edge( center, 1, core::kinematics::Edge::PEPTIDE );
    ft.add_edge( center, node1.second, core::kinematics::Edge::PEPTIDE );

    //Thinking about the loop closure problem
    //Need to skip 1 through center-1, add center through cutpoint
    //core::Size loop_count = 1; - use size of loop_vector_
    for (core::Size i = 1; i <= center-1; i++) {
        loop_for_residue_.push_back(0);
    }
    for (core::Size i = center; i < node1.second; i++) {
        loop_vector_.push_back(protocols::loops::Loop(i, node1.second+2, node1.second )); // this may be too aggressive
        loop_for_residue_.push_back(loop_vector_.size());
    }
    //skipping cutpoint for now
    loop_for_residue_.push_back(0);

    core::Size count = 1;
    //Creating the foldtree for given SS and previous loop
    for (size_t i = 2; i <= ss_stretch.size()-1; i++){
        //assuming loop size 3+ between SS elements
        std::cout << "starting " << i << std::endl;
        std::pair< core::Size, core::Size > node_i = ss_stretch[i];
        core::Size center_i = (node_i.first + node_i.second)/2;

        //std::pair< core::Size, core::Size > node_prev = input_vect[i-1];
        //core::Size loop_start = node_prev.second + 1; //prev node is SS, need loop
        //core::Size loop_end = node_i.first -1;
        //core::Size loop_center = (loop_start + loop_end)/2;

        //std::cout << "New jump " << center << " " << loop_center << " " << count << std::endl;
        //std::cout << "New edge " << loop_center << " " << loop_start << std::endl;
        //std::cout << "New edge " << loop_center << " " << loop_end << std::endl;
        //ft.add_edge( center, loop_center, count );
        //ft.add_edge( loop_center, loop_start, core::kinematics::Edge::PEPTIDE );
        //ft.add_edge( loop_center, loop_end, core::kinematics::Edge::PEPTIDE );

        //count++;

        std::cout << "New jump " << center << " " << center_i << " " << count << std::endl;
        std::cout << "New edge " << center_i << " " << node_i.first << std::endl;
        std::cout << "New edge " << center_i << " " << node_i.second << std::endl;
        ft.add_edge( center, center_i, count );
        ft.add_edge( center_i, node_i.first, core::kinematics::Edge::PEPTIDE );
        ft.add_edge( center_i, node_i.second, core::kinematics::Edge::PEPTIDE );

        count++;

        //ignoring actual cutpoints
        loop_for_residue_.push_back(0);
        for (core::Size j = node_i.first +1; j < center_i; j++) {
            loop_vector_.push_back(protocols::loops::Loop(j, node_i.first-2, node_i.first )); // this may be too aggressive
            loop_for_residue_.push_back(loop_vector_.size());
        }
        for (core::Size j = center_i; j < node_i.second; j++) {
            loop_vector_.push_back(protocols::loops::Loop(j, node_i.second+2, node_i.second )); // this may be too aggressive
            loop_for_residue_.push_back(loop_vector_.size());
        }
        loop_for_residue_.push_back(0);
    }

    //Making the end
    //assuming loop size 3+ between SS elements
    std::cout << "starting " << ss_stretch.size() << " (end)" << std::endl;
    std::pair< core::Size, core::Size > node_end = ss_stretch[ss_stretch.size()];
    core::Size center_end = (node_end.first + node_end.second)/2;

    std::cout << "New jump " << center << " " << center_end << " " << count << std::endl;
    std::cout << "New edge " << center_end << " " << node_end.first << std::endl;
    std::cout << "New edge " << center_end << " " << node_end.second << std::endl;
    ft.add_edge( center, center_end, count );
    ft.add_edge( center_end, node_end.first, core::kinematics::Edge::PEPTIDE );
    ft.add_edge( center_end, node_end.second, core::kinematics::Edge::PEPTIDE );
    //count++;

    //ignoring actual cutpoints
    loop_for_residue_.push_back(0);
    for (core::Size j = node_end.first +1; j < center_end; j++) {
        loop_vector_.push_back(protocols::loops::Loop(j, node_end.first-2, node_end.first )); // this may be too aggressive
        loop_for_residue_.push_back(loop_vector_.size());
    }
    for (core::Size j = center_end; j <= ss_stretch.size(); j++) {
        //loop_vector_.push_back(protocols::loops::Loop(j, node_i.second+2, node_i.second )); // this may be too aggressive
        loop_for_residue_.push_back(0);
    }
    //loop_for_residue_.push_back(0);

    ft_ = ft;

    //Making the end
//    std::pair< core::Size, core::Size > node_end = input_vect[input_vect.size()];
//
//    std::pair< core::Size, core::Size > node_prev = input_vect[input_vect.size()-1];
//    core::Size loop_start = node_prev.second + 1; //prev node is SS, need loop
//    core::Size loop_end = node_end.first -1;
//    core::Size loop_center = (loop_start + loop_end)/2;
//
//    std::cout << "New jump " << center << " " << loop_center << " " << count << std::endl;
//    std::cout << "New edge " << loop_center << " " << loop_start << std::endl;
//    std::cout << "New edge " << loop_center << " " << loop_end << std::endl;
//    ft.add_edge( center, loop_center, count );
//    ft.add_edge( loop_center, loop_start, core::kinematics::Edge::PEPTIDE );
//    ft.add_edge( loop_center, loop_end, core::kinematics::Edge::PEPTIDE );
//    count++;
//
//    core::Size center_end = (node_end.first + node_end.second)/2;
//    std::cout << "New jump " << center << " " << center_end << " " << count << std::endl;
//    std::cout << "New edge " << center_end << " " << node_end.first << std::endl;
//    std::cout << "New edge " << center_end << " " << input_string.size() << std::endl;
//    ft.add_edge( center, center_end, count );
//    ft.add_edge( center_end, node_end.first, core::kinematics::Edge::PEPTIDE );
//    //Need to get last element
//    ft.add_edge( center_end, input_string.size(), core::kinematics::Edge::PEPTIDE );

    //return ft;
}

//core::kinematics::FoldTree
void FoldTreeFromSS::fold_tree_from_ss(core::pose::Pose & pose) {
    core::scoring::dssp::Dssp dssp_struct = core::scoring::dssp::Dssp(pose);
    //core::scoring::dssp::Dssp dssp_mem = core::scoring::dssp::Dssp();
    std::string input_string = dssp_struct.get_dssp_secstruct();
    //return
    fold_tree_from_dssp_string(input_string);
}

}
}
