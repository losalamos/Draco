//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q1GaussLeg.hh
 * \author Kelly Thompson
 * \date   Wed Sep  1 09:35:03 2004
 * \brief  A class to encapsulate a 1D Gauss Legendre Quadrature set.
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_Q1GaussLeg_hh
#define quadrature_Q1GaussLeg_hh

#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class Q1GaussLeg
 * \brief A class to encapsulate a 1D Gauss Legendre Quadrature set.
 *
 * The client only needs to call QuadCreator::QuadCreate with the requested
 * quadrature set specified.  Since this quadrature set is inheireted from
 * the class Quadrature the client never needs to access this class directly.
 * The class overrides the virutal member functions contained in the class
 * Quadrature and contains data members that define the quadrature set.
 *
 * \sa Q1GaussLeg.cc for detailed descriptions.
 *
 * \example quadrature/test/tstQ1DGaussLeg.cc 
 *
 */
//===========================================================================//

class Q1DGaussLeg : public Quadrature
{
  public:

    // CREATORS

    // The default values for snOrder_ and norm_ were set in QuadCreator.
    Q1DGaussLeg( size_t snOrder_, double norm_, Quadrature::QIM qm_ );
    Q1DGaussLeg(); // disable default construction

    // ACCESSORS

    // These functions override the virtual member functions specifed in the
    // parent class Quadrature.
    
    size_t getNumOrdinates()   const { return numOrdinates; }
    void   display()        const;
    string name()           const { return "1D Gauss Legendre"; }
    string parse_name()     const { return "gauss legendre"; }
    size_t dimensionality() const { return 1; }
    size_t getSnOrder()     const { return snOrder; }
    Quadrature_Class getClass() const { return ONE_DIM; }

  private:

    // DATA
    size_t numOrdinates;  // == snOrder
};

} // end namespace rtt_quadrature

#endif // quadrature_Q1GaussLeg_hh

//---------------------------------------------------------------------------//
//              end of quadrature/Q1GaussLeg.hh
//---------------------------------------------------------------------------//
