//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/test/MatrixP1DiffUtils.hh
 * \author Randy M. Roberts
 * \date   Fri Jan 21 11:08:27 2000
 * \brief  Utilities to help test MatrixP1Diff
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMGDiffusionSolver_test_MatrixP1DiffUtils_hh__
#define __LAMGDiffusionSolver_test_MatrixP1DiffUtils_hh__

#include "../MatrixP1Diff.hh"
#include "LAMG/CompressedRowStorage.hh"

namespace rtt_LAMGDiffusionSolver_test
{
 
rtt_LAMG::CompressedRowStorage::DMat1
operator*(const rtt_LAMGDiffusionSolver::MatrixP1Diff &A,
	  const rtt_LAMG::CompressedRowStorage::DMat1 &x_in);

template<class FT>
void multiply(FT &b_out, const rtt_LAMGDiffusionSolver::MatrixP1Diff &A,
	      const FT &x_in)
{
    typedef rtt_LAMG::CompressedRowStorage::DMat1 DMat1;

    DMat1 x(x_in.size());
    std::copy(x_in.begin(), x_in.end(), x.begin());

    std::cout << "In templated multiply: x = " << std::endl;
    std::copy(x.begin(), x.end(),
	      std::ostream_iterator<DMat1::value_type>(std::cout, ","));
    std::cout << std::endl;

    DMat1 b = A * x;

    std::cout << "In templated multiply: b = " << std::endl;
    std::copy(b.begin(), b.end(),
	      std::ostream_iterator<DMat1::value_type>(std::cout, ","));
    std::cout << std::endl;
    
    std::copy(b.begin(), b.end(), b_out.begin());
}

template<class FT>
FT operator*(const rtt_LAMGDiffusionSolver::MatrixP1Diff &A, const FT &x_in)
{
    FT b(x_in.size());
    multiply(b, A, x_in);
    return b;
}
	      
} // end namespace rtt_LAMGDiffusionSolver_test

#endif                          // __LAMGDiffusionSolver_test_MatrixP1DiffUtils_hh__

//---------------------------------------------------------------------------//
//                              end of LAMGDiffusionSolver/test/MatrixP1DiffUtils.hh
//---------------------------------------------------------------------------//
