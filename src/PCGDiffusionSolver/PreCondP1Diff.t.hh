//----------------------------------*-C++-*----------------------------------//
// PreCondP1Diff.t.cc
// Randy M. Roberts
// Mon Sep 28 16:30:41 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "PreCondP1Diff.hh"
#include <stdexcept>
#include <string>
#include <sstream>

namespace rtt_PCGDiffusionSolver
{

 //---------------------------------------------------------------------------//
 // Solve QL * x = b for x.
 //---------------------------------------------------------------------------//

 template<class Matrix>
 void PreCondP1Diff<Matrix>::Left_PreCond(Mat1 &x, const Mat1 &b)
 {
     switch(method)
     {
     case 1:
	 // Here is Left Jacobi preconditioning ("diagonal scaling").
	 spMatrix->jacobiIteration(x, b);
	 break;
     case 0:
     case 2:
	 // This is Left "null preconditioning".
	 // x = b;
	 // break;
     default:
	 std::ostringstream ost;
	 ost << "Unrecognized Left preconditioning option: " << method;
	 throw std::runtime_error(ost.str());
     }
 }

 //---------------------------------------------------------------------------//
 // Solve QR * x = b for x.
 //---------------------------------------------------------------------------//

 template<class Matrix>
 void PreCondP1Diff<Matrix>::Right_PreCond(Mat1 &x, const Mat1 &b)
 {
     switch(method)
     {
     case 2:
	 // Here is Right Jacobi preconditioning ("diagonal scaling").
	 spMatrix->jacobiIteration(x, b);
	 break;
     case 0:
     case 1:
	 // This is Right "null preconditioning".
	 // x = b;
	 // break;
     default:
	 std::ostringstream ost;
	 ost << "Unrecognized Right preconditioning option: " << method;
	 throw std::runtime_error(ost.str());
     }
 }
 
} // end namespace rtt_PCGDiffusionSolver

//---------------------------------------------------------------------------//
//                              end of PreCondP1Diff.t.cc
//---------------------------------------------------------------------------//
