//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGrad/C4Reduction.hh
 * \author Shawn Pautz
 * \date   Mon May  7 14:57:07 2001
 * \brief  Provides C4 reduction operation for parallel calculations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGrad_C4Reduction_hh__
#define __ConjGrad_C4Reduction_hh__

#include "c4/global.hh"

namespace rtt_ConjGrad
{
 
//===========================================================================//
/*!
 * \class C4Reduction
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class C4Reduction 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
//     C4Reduction();
//     C4Reduction(const C4Reduction &rhs);
//     ~C4Reduction();

    // MANIPULATORS
    
//     C4Reduction& operator=(const C4Reduction &rhs);

    // ACCESSORS

    template <class T>
    T sum(const T &ret) const
    {
	T temp(ret);
	C4::gsum(temp);
	return temp;
    }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_ConjGrad

#endif                          // __ConjGrad_C4Reduction_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGrad/C4Reduction.hh
//---------------------------------------------------------------------------//
