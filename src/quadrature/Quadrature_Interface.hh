//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Quadrature_Interface.hh
 * \author Jae Chang
 * \date   Tue Jan 27 08:51:19 2004
 * \brief  Quadrature interface definitions
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
#ifndef quadrature_Quadrature_Interface_hh
#define quadrature_Quadrature_Interface_hh

// indirectly defines the CPP macro DLL_PUBLIC_quadrature
#include "ds++/Assert.hh"

extern "C" {

//===========================================================================//
/*!
 * \class quadrature_data
 * \brief Flattened quadrature data used for communicating with Fortran
 *        routines.
 */
//===========================================================================//

struct quadrature_data {
  int dimension;
  int type;
  int order;
  int azimuthal_order;
  int geometry;
  double *mu;
  double *eta;
  double *xi;
  double *weights;

  //! Default constructor for quadrature_data
  DLL_PUBLIC_quadrature quadrature_data()
      : dimension(0), type(0), order(0), azimuthal_order(0), geometry(0),
        mu(NULL), eta(NULL), xi(NULL), weights(NULL) { /* empty */
  }
};

//! An extern "C" interface to default constructor
DLL_PUBLIC_quadrature void init_quadrature(quadrature_data &quad);

//! Get quadrature data (eg. wts and cosines)
DLL_PUBLIC_quadrature void get_quadrature(quadrature_data &);

//! Ensure quadrature data is meaningful
DLL_PUBLIC_quadrature void check_quadrature_validity(const quadrature_data &);

} // end extern "C" block

#endif // quadrature_Quadrature_Interface_hh

//----------------------------------------------------------------------------//
// end of quadrature/Quadrature_Interface.hh
//----------------------------------------------------------------------------//
