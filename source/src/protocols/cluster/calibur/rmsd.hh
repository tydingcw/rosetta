// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file external/calibur/rmsd.hh
/// @author YK Ng & SC Li (kalngyk@gmail.com)

#ifndef external_calibur_rmsd_HH
#define external_calibur_rmsd_HH

#ifndef PYROSETTA
#endif
#if !defined(__WIN32__) && !defined(WIN32)
#endif
#include <vector>

#include <core/types.hh>

namespace protocols {
namespace cluster {
namespace calibur {

//#define __DEBUG_RMSD__
//#define __RMSD_SPEED_TEST__

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

#ifdef __RMSD_SPEED_TEST__
#ifndef PYROSETTA
#ifndef __WIN32__
inline core::Real
__timeval_difference(struct timeval * x, struct timeval * y) {
	core::Real elapsed;
	if (x->tv_usec < y->tv_usec) {
		int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
		y->tv_usec -= 1000000 * nsec;
		y->tv_sec += nsec;
	}
	if (x->tv_usec - y->tv_usec > 1000000) {
		int nsec = (x->tv_usec - y->tv_usec) / 1000000;
		y->tv_usec += 1000000 * nsec;
		y->tv_sec -= nsec;
	}

	elapsed = x->tv_sec - y->tv_sec;
	elapsed += (x->tv_usec - y->tv_usec)/(core::Real)1000000;
	return elapsed;
}

inline core::Real
_get_elapsed(int set_start) {
	static struct rusage last;
	struct rusage now;
	core::Real elapsed = 0;
	if (set_start) {
		getrusage(RUSAGE_SELF, &last);
	} else {
		getrusage(RUSAGE_SELF, &now);
		elapsed += __timeval_difference(&(now.ru_utime), &(last.ru_utime));
		elapsed += __timeval_difference(&(now.ru_stime), &(last.ru_stime));
		last = now;
	}
	return elapsed;
}
#endif
#endif
#endif

core::Real RMSD( std::vector<core::Real> & coords1, std::vector<core::Real> & coords2, int n, core::Real R[3][3] );

// Same as above, but skips computation of R
core::Real RMSD( std::vector<core::Real> & coords1, std::vector<core::Real> & coords2, int n );

void rotate( std::vector<core::Real> & coords, int n, core::Real R[3][3], core::Real * result );

core::Real fast_rmsd( std::vector<core::Real> & coords1, std::vector<core::Real> & coords2, int n);

}
}
}

#endif
