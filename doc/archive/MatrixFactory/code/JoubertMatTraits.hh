/*-----------------------------------*-C++-*---------------------------------*/
/* JoubertMatTraits.hh */
/* Randy M. Roberts */
/* Thu May 27 15:55:53 1999 */
/*---------------------------------------------------------------------------*/
/* @> */
/*---------------------------------------------------------------------------*/

#ifndef __MatrixFactory_JoubertMatTraits_hh__
#define __MatrixFactory_JoubertMatTraits_hh__

#include "JoubertMat.hh"
#include "MatrixFactoryTraits.hh"
#include "CRSMatrixRep.hh"

namespace rtt_MatrixFactory
{

// Specialize the MatrixFactoryTraits on a JoubertMat.

template<>
struct MatrixFactoryTraits<JoubertNS::JoubertMat>
{
    typedef JoubertNS::JoubertMat JoubertMat;

    // We have decided that any Draco-defined representation
    // that can be converted to a CRSMatrixRep can be used with
    // the "create" factory trait.  (This, of course, will require
    // the possible creation of a temporary CRSMatrixRep object.)
    // Anything that cannot be converted to a CRSMatrixRep will fail
    // at compile time.

    template<class T>
    static JoubertMat *create(const T &rep)
    {
	// The magic here is contained in the
	// "static_cast<CRSMatrixRep>(rep)" construct.
	//
	// There is another create utility overloaded to
	// take a const reference to a CRSMatrixRep object.
	// This overloaded create utility will be called by this
	// templated create utility.
	//
	// This means that a CRSMatrixRep object will be created from a
	// "T" object via type conversion.
	//
	// This conversion can be accomplished in the CRSMatrixRep
	// class via an explicit or non-explicit constructor taking
	// a const reference to a "T", or, alternately, in the "T" class
	// via an "operator CRSMatrixRep()" method.
	//
	// If this conversion is not defined, then this utility should
	// fail at compile time.  If it is defined then the "T"
	// representation is converted into a CRSMatrixRep representation.
	//
	// The CRSMatrixRep object is given to the overloaded create
	// utility of this traits class, whose job is to create the JoubertMat.

	return create(static_cast<CRSMatrixRep>(rep));
    }

    // Overload the "create" static method when given a CRSMatrixRep.
    
    static JoubertMat *create(const CRSMatrixRep &rep);
};

} // namespace rtt_MatrixFactory

#endif    /* __MatrixFactory_JoubertMatTraits_hh__ */

/*---------------------------------------------------------------------------*/
/*    end of JoubertMatTraits.hh */
/*---------------------------------------------------------------------------*/
