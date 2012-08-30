//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q3DLevelSym.hh
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \brief  A class to encapsulate a 3D Level Symmetric quadrature set.
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_Q3DLevelSym_hh
#define quadrature_Q3DLevelSym_hh

#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class Q3DLevelSym
 * \brief A class to encapsulate a 3D Level Symmetric quadrature set.
 *
 * The client only needs to call QuadCreator::QuadCreate with the requested
 * quadrature set specified.  Since this quadrature set is inheireted from
 * the class Quadrature the client never needs to access this class directly.
 * The class overrides the virutal member functions contained in the class
 * Quadrature and contains data members that define the quadrature set.
 *
 * \sa Q3DLevelSym.cc for detailed descriptions.
 * \example quadrature/test/tstQ3DLevelSym.cc 
 */
//===========================================================================//

class Q3DLevelSym : public Quadrature
{
  public:

    // CREATORS

    // The default values for snOrder_ and norm_ were set in QuadCreator.
    Q3DLevelSym( size_t snOrder_, double norm_, Quadrature::QIM qm_ );
    Q3DLevelSym();    // disable default construction

    // ACCESSORS

    // These functions override the virtual member functions specifed in the
    // parent class Quadrature.

    //! Returns the number of ordinates in the current quadrature set.
    size_t getNumOrdinates()   const { return numOrdinates; }
    //! Prints a short table containing the quadrature directions and weights.
    void display()       const;
    //! Returns the official name of the current quadrature set.
    string name()        const { return "3D Level Symmetric"; }
    //! Returns the input deck name of the current quadrature set.
    string parse_name()  const { return "level symmetric 3D"; }
    //! Returns the number of dimensions in the current quadrature set.
    size_t dimensionality() const { return 3; }
    //! Returns the order of the SN set.
    size_t getSnOrder()     const { return snOrder; }
    //! Returns the number of xi levels in the quadrature set.
    size_t getLevels()      const { return (snOrder+2)*snOrder/8; }
    //! Returns the class of quadrature.
    Quadrature_Class getClass() const { return THREE_DIM_TRIANGULAR; }

  private:

    // DATA
    size_t numOrdinates; // defaults to 24.
};

} // end namespace rtt_quadrature

#endif // quadrature_Q3DLevelSym_hh

//---------------------------------------------------------------------------//
//              end of quadrature/Q3DLevelSym.hh
//---------------------------------------------------------------------------//
