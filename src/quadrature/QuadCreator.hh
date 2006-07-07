//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadCreator.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:46:17 2000
 * \brief  Quadrature Creator class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __quadrature_QuadCreator_hh__
#define __quadrature_QuadCreator_hh__

#include "ds++/SP.hh"
#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class QuadCreator
 *
 * \brief A class to instantiate Quadrature objects.
 *
 * The generation of Quadrature objects can be implemented as a
 * "Parameterized Factory Method" (Booch).  This class acts as the virtual
 * creator for all quadrature objects.
 */

/*!
 * \example quadrature/test/tQuadrature.cc
 *
 * Example of Quadrature usage.  In this example the client must spcify the
 * quadrature type (must be one of the valid types specified by the
 * enumeration QuadCreator::Qid).  The SN order is optional and will be
 * defaulted to 4 is not specified.  The client may also specify a
 * normalization constant for the sum of of the direction weights.  If
 * unspecified, this normalization will default to 2, 2*PI or 4*PI for 1D, 2D 
 * or 3D quadratures, respectively.
 * 
 */
// revision history:
// -----------------
// 1.1) original
// 1.2) Implemented use of smart pointers (QuadCreator::QuadCreate now
//         returns a smartpointer instead of a normal pointer.)
//      Added/modified comments (both DOxygen and normal). 
//      Forced the default value for "norm" to be zero.  If it is zero
//         then "QuadCreate" will set norm to an appropriate default
//         based on the dimensionality of the quadrature set.
// 1.3) Renamed member function "QuadCreate" to start with a lower
//         case letter ("quadCreate").
// 
//===========================================================================//

class QuadCreator 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

  public:

    // DATA
    
    /*!
     * \brief A list of available quadrature types.
     *
     * Qid contains a list of identifiers that may be used to specify
     * the construction of a particular type of quadrature set (see
     * QuadCreate).  This list will grow as more types of quadrature are
     * added to this package. 
     */
    enum Qid { 
	SquareCL,    /*!< Square 2D Chebyshev-Legendre */
	GaussLeg,    /*!< 1D Gauss Legendre (arbitrary order). */
	Lobatto,     /*!< 1D Lobatto (arbitrary order). */
	DoubleGauss, /*!< 1D Double Gauss (arbitrary order). */
	LevelSym2D,  /*!< 2D Level Symmetric (even order between 2 and 24, inclusive). */
	LevelSym,    /*!< 3D Level Symmetric (even order between 2 and 24, inclusive). */
	Axial1D      /*!< 1D Axial used for filter sweeps */
    };

    // CREATORS

    // I'm not sure if this needs to be virtual or not.
    virtual rtt_dsxx::SP<Quadrature> quadCreate( Qid quad_type,
						 size_t sn_order = 4,
						 double norm = 0.0 );

};

} // end namespace rtt_quadrature

#endif // __quadrature_QuadCreator_hh__

//---------------------------------------------------------------------------//
//                              end of quadrature/QuadCreator.hh
//---------------------------------------------------------------------------//
