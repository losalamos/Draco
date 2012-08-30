//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q1DDoubleGauss.hh
 * \author James Warsa
 * \date   Fri Sep 16 15:45:26 2005
 * \brief  1D Double-Gauss Quadrature
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_Q1DDoubleGauss_hh
#define quadrature_Q1DDoubleGauss_hh

#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class Q1DDoubleGauss
 * \brief A class to encapsulate a 1D Lobatto Quadrature set.
 *
 * The client only needs to call QuadCreator::QuadCreate with the requested
 * quadrature set specified.  Since this quadrature set is inheireted from
 * the class Quadrature the client never needs to access this class directly.
 * The class overrides the virutal member functions contained in the class
 * Quadrature and contains data members that define the quadrature set.
 *
 * \sa Q1Lobatto.cc for detailed descriptions.
 *
 * \example quadrature/test/tstQ1DDoubleGauss.cc 
 *
 */
//===========================================================================//

class Q1DDoubleGauss : public Quadrature
{
  public:

    // CREATORS

    // The default values for snOrder_ and norm_ were set in QuadCreator.

    Q1DDoubleGauss( size_t snOrder_, double norm_, Quadrature::QIM qm_ );
    Q1DDoubleGauss(); // disable default construction

    // ACCESSORS

    // These functions override the virtual member functions specifed in the
    // parent class Quadrature.
    
    size_t getNumOrdinates()   const { return numOrdinates; }
    void   display()        const;
    string name()           const { return "1D Double-Gauss"; }
    string parse_name()     const { return "double gauss"; }
    size_t dimensionality() const { return 1; }
    size_t getSnOrder()     const { return snOrder; }
    Quadrature_Class getClass() const { return ONE_DIM; }

  private:

    // DATA
    size_t numOrdinates;  // == snOrder
};

} // end namespace rtt_quadrature

#endif // quadrature_Q1DDoubleGauss_hh

//---------------------------------------------------------------------------//
//              end of quadrature/Q1DDoubleGauss.hh
//---------------------------------------------------------------------------//
