//----------------------------------*-C++-*----------------------------------//
// P1Matrix.hh
// Randy M. Roberts
// Tue Jun  8 09:25:09 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __P1Diffusion_P1Matrix_hh__
#define __P1Diffusion_P1Matrix_hh__

#include "ds++/SP.hh"

namespace rtt_P1Diffusion
{
 
//===========================================================================//
// class P1Matrix - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
class P1Matrix 
{

    // NESTED CLASSES AND TYPEDEFS

     typedef typename MT::fcdsf fcdsf;
     typedef typename MT::ccsf ccsf;
     typedef typename MT::FieldConstructor FieldConstructor;

    // DATA
    
    const dsxx::SP<const ccsf> spADiagonal_m;
    const dsxx::SP<const fcdsf> spAOffDiagonal_m;
    const FieldConstructor &fCtor_m;
    
  public:

    // CREATORS
    
    P1Matrix(const FieldConstructor &fCtor_,
	     const dsxx::SP<const ccsf> &spADiag_,
	     const dsxx::SP<const fcdsf> &spAOffDiag_)
	: fCtor_m(fCtor_), spADiagonal_m(spADiag_),
	  spAOffDiagonal_m(spAOffDiag_)
    {
	// empty
    }
    
    // MANIPULATORS

    // ** none **
    
    // ACCESSORS

    const FieldConstructor &fCtor() const
    {
	return fCtor_m;
    }
    const dsxx::SP<const ccsf> spADiagonal() const
    {
	return spADiagonal_m;
    }
    const dsxx::SP<const fcdsf> spAOffDiagonal() const
    {
	return spAOffDiagonal_m;
    }
    
  private:
    
    // IMPLEMENTATION

    // ** none **
};

} // end namespace rtt_P1Diffusion

#endif                          // __P1Diffusion_P1Matrix_hh__

//---------------------------------------------------------------------------//
//                              end of P1Diffusion/P1Matrix.hh
//---------------------------------------------------------------------------//
