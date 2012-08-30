//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q2DTriChebyshevLegendre.hh
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

#ifndef quadrature_Q2DTriChebyshevLegendre_hh
#define quadrature_Q2DTriChebyshevLegendre_hh

#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class Q2DTriChebyshevLegendre
 * \brief
 *
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at http://www.doxygen.org.
 *
 * \sa Q2DTriChebyshevLegendre.cc for detailed descriptions.
 *
 * Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
 */
/*! 
 * \example quadrature/test/tstQ2DTriChebyshevLegendre.cc 
 * 
 * description of example
 */
//===========================================================================//

class Q2DTriChebyshevLegendre  : public Quadrature
{

  public:

    //! Disable default constructor.
    Q2DTriChebyshevLegendre( size_t snOrder_, double norm_, Quadrature::QIM qm_ );
    Q2DTriChebyshevLegendre();

    // ACCESSORS

    //! Returns the number of ordinates in the current quadrature set.
    size_t getNumOrdinates()   const { return numOrdinates; }
    //! Prints a short table containing the quadrature directions and weights.
    void display()       const;
    //! Returns the official name of the current quadrature set.
    string name()        const { return "2D Tri Chebyshev Legendre"; }
    //! Returns the input deck name of the current quadrature set.
    string parse_name()  const { return "tri cl"; }
    //! Returns the number of dimensions in the current quadrature set.
    size_t dimensionality() const { return 2; }
    //! Returns the order of the SN set.
    size_t getSnOrder()     const { return snOrder; }
    //! Returns the number of eta levels in the quadrature set.
    size_t getLevels()      const { return snOrder; }
    //! Returns the class of quadrature.
    Quadrature_Class getClass() const { return TWO_DIM_TRIANGULAR; }

    // MANIPULATORS
    void sortOrdinates(void);

  private:

    size_t numOrdinates;
};

} // end namespace rtt_quadrature

#endif // quadrature_Q2DTriChebyshevLegendre_hh

//---------------------------------------------------------------------------//
//              end of quadrature/Q2DTriChebyshevLegendre.hh
//---------------------------------------------------------------------------//
