/*-----------------------------------*-C-*-----------------------------------*/
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

template<>
struct MatrixFactoryTraits<JoubertMat::JoubertMat>
{
    typedef JoubertMat::JoubertMat Matrix;

    template<class T>
    static Matrix *create(const T &rep)
    {
	// The magic here is contained in the
	// "static_cast<CRSMatrixRep>(rep)" construct.
	//
	// There is create utility takes a const reference to
	// a CRSMatrixRep object.
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
	// utility of this traits class, whose job is to create the Matrix.

	return create(static_cast<CRSMatrixRep>(rep));
    }
    
    static Matrix *create(const CRSMatrixRep &rep);
};

} // namespace rtt_MatrixFactory

#endif    /* __MatrixFactory_JoubertMatTraits_hh__ */

/*---------------------------------------------------------------------------*/
/*    end of JoubertMatTraits.hh */
/*---------------------------------------------------------------------------*/
