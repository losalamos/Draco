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
	GaussLeg,    /*!< 1D Gauss Legendre (arbitrary order). */
	LevelSym2D,  /*!< 2D Level Symmetric (even order between 2 and 24, inclusive). */
	LevelSym     /*!< 3D Level Symmetric (even order between 2 and 24, inclusive). */
    };

    // CREATORS

    /*!
     * \brief quadCreate constructs a Quadrature object.
     *
     * The Quad creator only requires 1 parameter -- the quadrature
     * identifier (see QuadCreator::Qid).  The two addtional parameters can
     * optionally be used to specify the sn_order and a normalization for the 
     * quadrature weights.  The sn_order defaults to 4 and the default value
     * for normalization constant varies with the dimensionality of the
     * quadrature set (2, 2*pi or 4*pi for 1D, 2D or 3D sets).
     *
     * Another parameter may need to be added to this constructor to specify
     * the number of dimensions requested.  Currently Qid directly specifies
     * the dimensionality of the quadrature set.
     *
     * \param quad_type An identifier that specifies the type of quadrature
     *                  to construct.
     * \param sn_order  The SN order for the constructed quadrature
     *                  set. (Default: 4)
     * \param norm      The sum of the quadrature weights are forced to sum
     *                  to this value. (Default: 2, 2*pi or 4*pi based on the 
     *                  dimensionality of the quadrature set.)
     * \return Smart pointer to a quadrature object.
     */
    // I'm not sure if this needs to be virtual or not.
    // norm and sn_order are optional variables.  sn_order always defaults to 
    // 4 but the default value for norm is based on the dimensionality of the 
    // quadrature set and its default value is set in the member function
    virtual rtt_dsxx::SP<Quadrature> 
          quadCreate( Qid quad_type, int sn_order = 4, double norm = 0.0 );
    //    QuadCreator(const QuadCreator &rhs);
    //   ~QuadCreator();

    // MANIPULATORS
    
    //    QuadCreator& operator=(const QuadCreator &rhs);

    // ACCESSORS

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_quadrature

#endif                          // __quadrature_QuadCreator_hh__

//---------------------------------------------------------------------------//
//                              end of quadrature/QuadCreator.hh
//---------------------------------------------------------------------------//
