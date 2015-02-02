//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Quadrature_Interface.hh
 * \author Jae Chang
 * \date   Tue Jan 27 08:51:19 2004
 * \brief  Quadrature interface definitions
 * \note   Copyright (C) 2004-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: Arguments.hh 7090 2015-01-07 17:01:48Z kellyt $
//---------------------------------------------------------------------------//

#ifndef quadrature_Quadrature_Interface_hh
#define quadrature_Quadrature_Interface_hh

extern "C"
{

struct quadrature_data
{
    int dimension;
    int type;
    int order;
    int azimuthal_order;
    int geometry;
    double *mu;
    double *eta;
    double *xi;
    double *weights;
    quadrature_data();

};

    //! An extern "C" interface to default constructor
    void init_quadrature(quadrature_data&);

    //! Get quadrature data (eg. wts and cosines)
    void get_quadrature(quadrature_data&);

}
#endif //quadrature_Quadrature_Interface_hh
