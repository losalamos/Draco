//----------------------------------*-C++-*----------------------------------//
// MatrixFactoryTraits.hh
// Randy M. Roberts
// Tue Jun  8 09:18:38 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __traits_MatrixFactoryTraits_hh__
#define __traits_MatrixFactoryTraits_hh__

namespace rtt_traits
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

    // DATA
    
  public:

    // STATIC CLASS METHODS

    template<class T>
    static Matrix *create(const T &rep)
    {
	// You should be specializing this class.
	// BogusMethod is being used to trigger a compilation
	// error.
	
	return MatrixFactoryTraits<Matrix>::BogusMethod(rep);
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

#endif                          // __traits_MatrixFactoryTraits_hh__

//---------------------------------------------------------------------------//
//                              end of traits/MatrixFactoryTraits.hh
//---------------------------------------------------------------------------//
