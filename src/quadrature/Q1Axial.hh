
//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q1Axial.hh
 * \author Jae H Chang
 * \date   Thu Oct  7 14:30:08 2004
 * \brief  A class to encapsulate a 1D Axial quadrature set.  Used in
 *         stretched and filter TSA.
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_Q1Axial_hh
#define quadrature_Q1Axial_hh

#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class Q1Axial
 * \brief A class to encapsulate a 1D Axial quadrature set.
 *
 * The client only needs to call QuadCreator::QuadCreate with the requested
 * quadrature set specified.  Since this quadrature set is inheireted from
 * the class Quadrature the client never needs to access this class directly.
 * The class overrides the virutal member functions contained in the class
 * Quadrature and contains data members that define the quadrature set.
 *
 * \sa Q1Axial.cc for detailed descriptions.
 *
 * \example quadrature/test/tstQ1Axial.cc 
 * 
 * 
 */
//===========================================================================//

class Q1Axial : public Quadrature
{
  public:

     // CREATORS
    
    //! constructors.
    Q1Axial( size_t snOrder_, double norm_, Quadrature::QIM qm_ );
    Q1Axial(); // disable default consturtion

     // ACCESSORS

    // These functions override the virtual member functions specifed in the
    // parent class Quadrature.
    
    size_t getNumOrdinates()   const { return numOrdinates; }
    void   display()        const;
    string name()           const { return "1D Axial"; }
    string parse_name()     const { return "axial"; }
    size_t dimensionality() const { return 1; }
    size_t getSnOrder()     const { return snOrder; }
    Quadrature_Class getClass() const { return ONE_DIM; }

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA
    size_t numOrdinates;  // == snOrder
};

} // end namespace rtt_quadrature

#endif // quadrature_Q1Axial_hh

//---------------------------------------------------------------------------//
//              end of quadrature/Q1Axial.hh
//---------------------------------------------------------------------------//
