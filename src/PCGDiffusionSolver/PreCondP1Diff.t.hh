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

namespace rtt_PCGDiffusionSolver
{

 //---------------------------------------------------------------------------//
 // Solve QL * x = b for x.
 //---------------------------------------------------------------------------//

 template<class Matrix>
 void PreCondP1Diff<Matrix>::Left_PreCond(Mat1 &x, Mat1 &b)
 {
     switch(method)
     {
     case 0:
	 // This is "null preconditioning".
	 x = b;
	 break;
     case 1:
	 // Here is Jacobi preconditioning ("diagonal scaling").
	 spMatrix->jacobiIteration(x, b);
	 break;
     default:
	 throw std::runtime_error("Unrecognized preconditioning option.");
     }
 }

 //---------------------------------------------------------------------------//
 // Solve QR * x = b for x.
 //---------------------------------------------------------------------------//

 template<class Matrix>
 void PreCondP1Diff<Matrix>::Right_PreCond(Mat1 &x, Mat1 &b)
 {
    x = b;
 }
 
} // end namespace rtt_PCGDiffusionSolver

//---------------------------------------------------------------------------//
//                              end of PreCondP1Diff.t.cc
//---------------------------------------------------------------------------//
