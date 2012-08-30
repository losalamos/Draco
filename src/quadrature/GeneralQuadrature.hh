//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/GeneralQuadrature.hh
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

#ifndef quadrature_GeneralQuadrature_hh
#define quadrature_GeneralQuadrature_hh

#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class GeneralQuadrature
 * \brief A class to encapsulate a generic user supplied quadrature set.
 *
 * \sa GeneralQuadrature.cc for detailed descriptions.
 * \example quadrature/test/tstGeneralQuadrature.cc 
 */
//===========================================================================//

class GeneralQuadrature : public Quadrature
{
  public:

    // CREATORS

    // The default values for snOrder_ and norm_ were set in QuadCreator.
    GeneralQuadrature( size_t snOrder_,
                       double norm_,
                       Quadrature::QIM qm,
                       std::vector<double>  const & mu_,
                       std::vector<double>  const & eta_,
                       std::vector<double>  const & xi_,
                       std::vector<double>  const & wt_,
                       size_t levels,
                       size_t dim,
                       std::string  const & quadratureName,
                       Quadrature_Class quadratureClass);

    GeneralQuadrature();    // disable default construction

    // ACCESSORS

    // These functions override the virtual member functions specifed in the
    // parent class Quadrature.

    //! Returns the number of ordinates in the current quadrature set.
    size_t getNumOrdinates()   const { return numOrdinates; }
    //! Prints a short table containing the quadrature directions and weights.
    void display()       const;
    //! Returns the official name of the current quadrature set.
    string name()        const { return quadratureName; }
    //! Returns the official name of the current quadrature set.
    string parse_name()  const { return quadratureName; }
    //! Returns the number of dimensions in the current quadrature set.
    size_t dimensionality() const { return numDims; }
    //! Returns the order of the SN set.
    size_t getSnOrder()     const { return snOrder; }
    //! Returns the number of xi levels in the quadrature set.
    size_t getLevels()      const { return numLevels; }

    Quadrature_Class getClass() const {return quadratureClass; };

  private:

    // DATA
    size_t const numOrdinates; // defaults to 24.
    size_t const numLevels;
    size_t const numDims;
    std::string const quadratureName;
    Quadrature_Class quadratureClass; 
};

} // end namespace rtt_quadrature

#endif // quadrature_GeneralQuadrature_hh

//---------------------------------------------------------------------------//
//              end of quadrature/GeneralQuadrature.hh
//---------------------------------------------------------------------------//
