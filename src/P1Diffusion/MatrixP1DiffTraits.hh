//----------------------------------*-C++-*----------------------------------//
// MatrixP1DiffTraits.hh
// Randy M. Roberts
// Tue Jun  8 09:50:15 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __P1Diffusion_MatrixP1DiffTraits_hh__
#define __P1Diffusion_MatrixP1DiffTraits_hh__

#include "traits/MatrixFactoryTraits.hh"
#include "MatrixP1Diff.hh"
#include "P1Matrix.hh"

namespace rtt_traits
{
 
//===========================================================================//
// class MatrixP1DiffTraits - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
class MatrixFactoryTraits<rtt_P1Diffusion::MatrixP1Diff<MT> >
{

    // NESTED CLASSES AND TYPEDEFS

    typedef rtt_P1Diffusion::MatrixP1Diff<MT> Matrix;

    // DATA
    
  public:

    // STATIC CLASS METHODS

    static Matrix *create(const rtt_P1Diffusion::P1Matrix<MT> &rep)
    {
	return new Matrix(rep.fCtor(), rep.spADiagonal(),
			  rep.spAOffDiagonal());
    }

    // CREATORS
    
    // ** none **

    // MANIPULATORS
    
    // ** none **

    // ACCESSORS

  private:
    
    // IMPLEMENTATION

    // ** none **
};

} // end namespace rtt_traits

#endif                          // __P1Diffusion_MatrixP1DiffTraits_hh__

//---------------------------------------------------------------------------//
//                              end of P1Diffusion/MatrixP1DiffTraits.hh
//---------------------------------------------------------------------------//
