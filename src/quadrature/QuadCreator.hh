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

#include <map>

#include "ds++/SP.hh"
#include "parser/Token_Stream.hh"
#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class QuadCreator
 *
 * \brief A class to instantiate Quadrature objects.
 *
 * The generation of Quadrature objects can be implemented as a "Parameterized
 * Factory Method" (Booch).  This class acts as the virtual creator for all
 * quadrature objects.
 */

/*!
 * \example quadrature/test/tQuadrature.cc
 *
 * Example of Quadrature usage.  In this example the client must spcify the
 * quadrature type (must be one of the valid types specified by the
 * enumeration QuadCreator::Qid).  The SN order is optional and will be
 * defaulted to 4 is not specified.  The client may also specify a
 * normalization constant for the sum of of the direction weights.  If
 * unspecified, this normalization will default to 1 (so that the zeroth
 * moment is equal to the mean over the sphere.)
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
// 1.17) Provide creator function for rtt_parser::Token_Stream.
// 
//===========================================================================//

class QuadCreator 
{
  public:

    // ENUMERATIONS
    
    /*!
     * \brief A list of available quadrature types.
     *
     * Qid contains a list of identifiers that may be used to specify
     * the construction of a particular type of quadrature set (see
     * QuadCreate).  This list will grow as more types of quadrature are
     * added to this package. 
     */
    enum Qid { 
	TriCL2D,     /*!< Triangular 2D Chebyshev-Legendre */
	SquareCL,    /*!< Square 2D Chebyshev-Legendre */
	GaussLeg,    /*!< 1D Gauss Legendre (arbitrary order). */
	Lobatto,     /*!< 1D Lobatto (arbitrary order). */
	DoubleGauss, /*!< 1D Double Gauss (arbitrary order). */
	LevelSym2D,  /*!< 2D Level Symmetric (even order between 2 and 24, inclusive). */
	LevelSym,    /*!< 3D Level Symmetric (even order between 2 and 24, inclusive). */
	TriCL,       /*!< Triangular 3D Chebyshev-Legendre */
	Axial1D      /*!< 1D Axial used for filter sweeps */
    };

    // TYPEDEFS

    typedef std::map< std::string, Qid > qidm;
    typedef std::map< Qid, double > normmap;

//    rtt_mesh_element::Geometry parsed_geometry;

    QuadCreator(void) : Qid_map( createQidMap() )
    { /* empty */ }

    virtual ~QuadCreator(void) { /* empty */ }
    
    // CREATORS

    // I'm not sure if this needs to be virtual or not.
    virtual rtt_dsxx::SP<Quadrature> quadCreate( Qid quad_type,
						 size_t sn_order = 4,
						 double norm = 0.0,
                                                 Quadrature::QIM interpModel = Quadrature::SN );

    rtt_dsxx::SP<Quadrature> quadCreate( rtt_parser::Token_Stream &tokens );

  private:

    // Functions
    qidm createQidMap(void) const;
        
    // DATA
    qidm    const Qid_map;    
};

} // end namespace rtt_quadrature

#endif // __quadrature_QuadCreator_hh__

//---------------------------------------------------------------------------//
//                              end of quadrature/QuadCreator.hh
//---------------------------------------------------------------------------//
