//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGrad/DefaultScalarMulties.hh
 * \author Randy M. Roberts
 * \date   Mon Apr 24 09:05:41 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGrad_DefaultScalarMulties_hh__
#define __ConjGrad_DefaultScalarMulties_hh__

#include <functional>

namespace rtt_ConjGrad
{
 
//===========================================================================//
/*!
 * \class DefaultScalarMulties
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class Arg1, class Arg2, class Result=Arg2>
class DefaultScalarMultiplies : public std::binary_function<Arg1,Arg2,Result>
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    //Defaulted: DefaultScalarMultiplies();
    //Defaulted:  DefaultScalarMultiplies(const DefaultScalarMultiplies &rhs);
    //Defaulted: ~DefaultScalarMultiplies();

    // MANIPULATORS
    
    //Defaulted: DefaultScalarMultiplies& operator=(
    //Defaulted:      const DefaultScalarMultiplies &rhs);

    // ACCESSORS

    Result operator()(Arg1 arg1, Arg2 arg2) const
    {
	return arg1*arg2;
    }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_ConjGrad

#endif                          // __ConjGrad_DefaultScalarMultiplies_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGrad/DefaultScalarMultiplies.hh
//---------------------------------------------------------------------------//
