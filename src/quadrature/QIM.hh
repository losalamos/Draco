//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/QIM.hh
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Definition of QIM enumeration
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------------------//
// $Id: QIM.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#ifndef quadrature_QIM_hh
#define quadrature_QIM_hh

namespace rtt_quadrature
{

//! Quadrature Interpolation Model: specifies how to compute the
//! Discrete-to-Moment operator. This is used by input parsers.

enum QIM 
{
    SN,  /*!< Use the standard SN method. */
    GQ,  /*!< Use Morel's Galerkin Quadrature method. */
    SVD  /*!< Let M be an approximate inverse of D. */
};

} // end namespace rtt_quadrature

#endif // quadrature_QIM_hh

//---------------------------------------------------------------------------------------//
//              end of quadrature/QIM.hh
//---------------------------------------------------------------------------------------//
