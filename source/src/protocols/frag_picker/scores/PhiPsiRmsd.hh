// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/PhiPsiRmsd.hh
/// @brief  Object that scores a fragment by its crmsd to the native
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_PhiPsiRmsd_hh
#define INCLUDED_protocols_frag_picker_scores_PhiPsiRmsd_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/PhiPsiTalosIO.fwd.hh>

#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// mini
#include <core/types.hh>

#include <ObjexxFCL/FArray1D.hh>

#include <utility/vector1.hh>

#include <core/pose/Pose.fwd.hh> // AUTO IWYU For PoseOP


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<core::Real> > Matrix;

/// @brief  scores a fragment by the root mean square deviation of Phi and Psi angles.
class PhiPsiRmsd: public CachingScoringMethod {
public:

	/// @brief  creates a Phi-Psi-based scoring function.
	/// @details Phi-Psi angles from a fragment will be compared to relevant angles in a given pose, which should have the same number of residues a the query sequence
	PhiPsiRmsd(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,
		core::pose::PoseOP reference_pose);
	PhiPsiRmsd(core::Size, core::Real, bool, PhiPsiTalosIO&);
	void do_caching(VallChunkOP) override;
	void clean_up() override;
	bool score(FragmentCandidateOP, FragmentScoreMapOP) override;
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP) override;

	/// @brief  Print debugging informations about a score for a given fragment
	bool describe_score(FragmentCandidateOP, FragmentScoreMapOP, std::ostream&);

private:
	std::string cached_scores_id_;
	core::Size n_atoms_;
	ObjexxFCL::FArray1D_double chunk_phi_;
	ObjexxFCL::FArray1D_double chunk_psi_;
	ObjexxFCL::FArray1D_double query_phi_;
	ObjexxFCL::FArray1D_double query_psi_;
	utility::vector1<bool> existing_data_;
	void read_talos_phi_psi(std::string const&);
};

/// @brief  Matker class that produces a new PhiPsiRmsd object
class MakePhiPsiRmsd: public MakeFragmentScoringMethod {
public:

	MakePhiPsiRmsd() :
		MakeFragmentScoringMethod("PhiPsiRmsd") {
	}

	FragmentScoringMethodOP make(core::Size, core::Real, bool, FragmentPickerOP, std::string const &) override;
};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_PhiPsiRmsd_HH
