//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   quadrature/Quadrature_Interface.hh
 * \author Jae Chang
 * \date   Tue Jan 27 08:51:19 2004
 * \brief  Quadrature interface definitions
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//
#ifndef rtt_quadrature_Quadrature_Interface_hh
#define rtt_quadrature_Quadrature_Interface_hh

extern "C" {

//============================================================================//
/*!
 * \class quadrature_data
 * \brief Flattened quadrature data used for communicating with Fortran
 *        routines.
 */
//============================================================================//

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
  quadrature_data()
      : dimension(0), type(0), order(0), azimuthal_order(0), geometry(0),
        mu(nullptr), eta(nullptr), xi(nullptr), weights(nullptr) { /* empty */
  }
};

//! An extern "C" interface to default constructor
void init_quadrature(quadrature_data &quad);

//! Get quadrature data (e.g.: weights and cosines)
void get_quadrature(quadrature_data &);

//! Ensure quadrature data is meaningful
void check_quadrature_validity(const quadrature_data &);

} // end extern "C" block

#endif // quadrature_Quadrature_Interface_hh

//----------------------------------------------------------------------------//
// end of quadrature/Quadrature_Interface.hh
//----------------------------------------------------------------------------//
