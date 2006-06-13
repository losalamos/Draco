//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q2DSquareChebyshevLegendre.hh
 * \author James Warsa
 * \date   Fri Jun  9 13:52:25 2006
 * \brief  
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_Q2DSquareChebyshevLegendre_hh
#define quadrature_Q2DSquareChebyshevLegendre_hh

#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class Q2DSquareChebyshevLegendre
 * \brief
 *
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at http://www.doxygen.org.
 *
 * \sa Q2DSquareChebyshevLegendre.cc for detailed descriptions.
 *
 * Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
 */
/*! 
 * \example quadrature/test/tstQ2DSquareChebyshevLegendre.cc 
 * 
 * description of example
 */
//===========================================================================//

class Q2DSquareChebyshevLegendre  : public Quadrature
{

  public:

    //! Disable default constructor.
    Q2DSquareChebyshevLegendre( size_t snOrder_, double norm_ );
    Q2DSquareChebyshevLegendre();

    // ACCESSORS

    //! Returns the number of angles in the current quadrature set.
    size_t getNumAngles()   const { return numAngles; }
    //! Prints a short table containing the quadrature directions and weights.
    void display()       const;
    //! Returns the official name of the current quadrature set.
    string name()        const { return "2D Square Chebyshev Legendre"; }
    //! Returns the number of dimensions in the current quadrature set.
    size_t dimensionality() const { return 2; }
    //! Returns the order of the SN set.
    size_t getSnOrder()     const { return snOrder; }
    //! Returns the number of eta levels in the quadrature set.
    size_t getLevels()      const { return snOrder; }

  private:

    size_t numAngles;
};

} // end namespace rtt_quadrature

#endif // quadrature_Q2DSquareChebyshevLegendre_hh

//---------------------------------------------------------------------------//
//              end of quadrature/Q2DSquareChebyshevLegendre.hh
//---------------------------------------------------------------------------//
