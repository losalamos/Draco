//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGrad/SerialReduction.hh
 * \author Shawn Pautz
 * \date   Mon May  7 14:43:46 2001
 * \brief  Provides reduction operation for serial calculations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGrad_SerialReduction_hh__
#define __ConjGrad_SerialReduction_hh__

namespace rtt_ConjGrad
{
 
//===========================================================================//
/*!
 * \class SerialReduction
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class SerialReduction 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
//     SerialReduction();
//     SerialReduction(const SerialReduction &rhs);
//     ~SerialReduction();

    // MANIPULATORS
    
//     SerialReduction& operator=(const SerialReduction &rhs);

    // ACCESSORS

    template <class T>
    T sum(const T &ret) const
    {
	return ret;
    }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_ConjGrad

#endif                          // __ConjGrad_SerialReduction_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGrad/SerialReduction.hh
//---------------------------------------------------------------------------//
