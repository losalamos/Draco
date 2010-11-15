//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/gaussj.hh
 * \author Kent Budge
 * \brief  Solve a system of equations by Gaussian elimination.
 * \note   © Copyright 2006-2007 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef linear_gaussj_hh
#define linear_gaussj_hh

namespace rtt_linear
{

//! Solve a system of linear equations.  Single-subscript computed-index
//! version. 
template<class RandomContainer>
void gaussj(RandomContainer &A,
            unsigned n,
	    RandomContainer &b,
	    unsigned m);

//! Solve a system of linear equations. Double-subscript version.
template<class DoubleRandomContainer, class RandomContainer>
void gaussj(DoubleRandomContainer &A,
	    RandomContainer &b);

} // end namespace rtt_linear

#endif // linear_gaussj_hh

//---------------------------------------------------------------------------//
//              end of linear/gaussj.hh
//---------------------------------------------------------------------------//
