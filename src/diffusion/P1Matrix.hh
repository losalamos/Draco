//----------------------------------*-C++-*----------------------------------//
// P1Matrix.hh
// Randy M. Roberts
// Tue Jun  8 09:25:09 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __diffusion_P1Matrix_hh__
#define __diffusion_P1Matrix_hh__

#include "ds++/SP.hh"

namespace rtt_diffusion
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
    
    ccsf ADiagonal_m;
    fcdsf AOffDiagonal_m;
    const FieldConstructor &fCtor_m;
    
  public:

    // CREATORS
    
    P1Matrix(const FieldConstructor &fCtor_, const ccsf &ADiag_,
	     const fcdsf &AOffDiag_)
	: fCtor_m(fCtor_), ADiagonal_m(ADiag_),
	  AOffDiagonal_m(AOffDiag_)
    {
	// Empty
    }
    
    // MANIPULATORS

    void jacobiScale()
    {
	fcdsf DiagOnFaces(fCtor());
	MT::gather(DiagOnFaces, ADiagonal_m, MT::OpAssign());

	fcdsf DiagAcrossFaces(fCtor());
	MT::swap_faces(DiagAcrossFaces, DiagOnFaces, 1.0);

	AOffDiagonal_m = AOffDiagonal_m / sqrt(fabs(DiagAcrossFaces
						    * DiagOnFaces));
	ADiagonal_m = ADiagonal_m/ fabs(ADiagonal_m);
    }
 
    
    // ACCESSORS

    const FieldConstructor &fCtor() const
    {
	return fCtor_m;
    }
    const ccsf &diagonal() const
    {
	return ADiagonal_m;
    }
    const fcdsf &offDiagonal() const
    {
	return AOffDiagonal_m;
    }
    
  private:
    
    // IMPLEMENTATION

    // ** none **
};

} // end namespace rtt_diffusion

#endif                          // __diffusion_P1Matrix_hh__

//---------------------------------------------------------------------------//
//                              end of diffusion/P1Matrix.hh
//---------------------------------------------------------------------------//
