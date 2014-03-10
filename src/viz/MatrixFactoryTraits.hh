//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   viz/MatrixFactoryTraits.hh
 * \author Randy M. Roberts
 * \date   Tue Jun  8 09:18:38 1999
 * \brief  
 * \note   Copyright (C) 1998-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------------------//

#ifndef __viz_MatrixFactoryTraits_hh__
#define __viz_MatrixFactoryTraits_hh__

namespace rtt_viz
{
 
//===========================================================================//
// class MatrixFactoryTraits - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class Matrix>
class MatrixFactoryTraits 
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    // This is a bogus typedef.
    // Create your own, as needed.
    
    typedef int PreComputedState;

  private:

    // DATA
    
  public:

    // STATIC CLASS METHODS

    template<class T>
    static PreComputedState preComputeState(const T &rep)
    {
	// You should be specializing this class.
	// BogusMethod is being used to trigger a compilation
	// error.
	
	return MatrixFactoryTraits<Matrix>::BogusMethod1(rep);
    }

    template<class T>
    static Matrix *create(const T &rep, const PreComputedState &state)
    {
	// You should be specializing this class.
	// BogusMethod is being used to trigger a compilation
	// error.
	
	return MatrixFactoryTraits<Matrix>::BogusMethod2(rep, state);
    }

    // CREATORS
    
    // ** none **

    // MANIPULATORS
    
    // ** none **

    // ACCESSORS

    // ** none **

  private:
    
    // IMPLEMENTATION

    // ** none **
};

} // end namespace rtt_traits

#endif // __viz_MatrixFactoryTraits_hh__

//---------------------------------------------------------------------------//
// end of viz/MatrixFactoryTraits.hh
//---------------------------------------------------------------------------//
