//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGrad/DefaultConvCrit.hh
 * \author Randy M. Roberts
 * \date   Mon Apr 24 11:06:42 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGrad_DefaultConvCrit_hh__
#define __ConjGrad_DefaultConvCrit_hh__

#include <iostream>

namespace rtt_ConjGrad
{
 
//===========================================================================//
/*!
 * \class DefaultConvCrit
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class Field, class Norm>
class DefaultConvCrit 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename Field::value_type value_type;

    // DATA

    Norm norm;
    int maxIters;
    value_type eps;
    mutable int curIter;
    
  public:

    // CREATORS
    
    DefaultConvCrit(Norm norm_in, int maxIters_in, value_type eps_in)
	: norm(norm_in), maxIters(maxIters_in), eps(eps_in),
	  curIter(0)
    {
	/* empty */
    }
    
    //Defaulted: DefaultConvCrit(const DefaultConvCrit &rhs);
    //Defaulted: ~DefaultConvCrit();

    // MANIPULATORS
    
    //Defaulted: DefaultConvCrit& operator=(const DefaultConvCrit &rhs);

    // ACCESSORS

    bool operator()(int iter, const Field &r, const Field &x,
		    const Field &b) const
    {
	curIter = iter;
	if (iter >= maxIters)
	{
	    std::cerr << "WARNING: "
		      << " Maximum number of iterations, "
		      << maxIters << ", exceeded."
		      << std::endl;
	    return true;
	}

	return norm(r) <= eps * norm(b);
    }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_ConjGrad

#endif                          // __ConjGrad_DefaultConvCrit_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGrad/DefaultConvCrit.hh
//---------------------------------------------------------------------------//
