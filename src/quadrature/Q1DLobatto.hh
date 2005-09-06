//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /Q1DLobatto.hh
 * \author James Warsa
 * \date   Fri Sep  2 10:30:02 2005
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef _Q1DLobatto_hh
#define _Q1DLobatto_hh

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class Q1Lobatto
 * \brief A class to encapsulate a 1D Gauss Legendre Quadrature set.
 *
 * The client only needs to call QuadCreator::QuadCreate with the requested
 * quadrature set specified.  Since this quadrature set is inheireted from
 * the class Quadrature the client never needs to access this class directly.
 * The class overrides the virutal member functions contained in the class
 * Quadrature and contains data members that define the quadrature set.
 *
 * \sa Q1Lobatto.cc for detailed descriptions.
 *
 * \example quadrature/test/tstQ1DLobatto.cc 
 *
 */
//===========================================================================//

class Q1DLobatto : public Quadrature
{
  public:

    // CREATORS

    // The default values for snOrder_ and norm_ were set in QuadCreator.
    Q1DLobatto( size_t snOrder_, double norm_ );
    Q1DLobatto(); // disable default construction

    // ACCESSORS

    // These functions override the virtual member functions specifed in the
    // parent class Quadrature.
    
    size_t getNumAngles()   const { return numAngles; }
    void   display()        const;
    string name()           const { return "1D Lobatto"; }
    size_t dimensionality() const { return 1; }
    size_t getSnOrder()     const { return snOrder; }

  private:

    // DATA
    size_t numAngles;  // == snOrder
};

} // end namespace rtt_quadrature

#endif // _Q1DLobatto_hh

//---------------------------------------------------------------------------//
//              end of /Q1DLobatto.hh
//---------------------------------------------------------------------------//
