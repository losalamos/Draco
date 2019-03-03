//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh_element/Geometry.hh
 * \author Kent Budge
 * \date   Tue Dec 21 14:28:56 2004
 * \brief
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef mesh_element_Geometry_hh
#define mesh_element_Geometry_hh

namespace rtt_mesh_element {

//! Enumerates supported geometries.
enum Geometry {
  AXISYMMETRIC, //!< 2D (cylindrical) R-Z
  SPHERICAL,    //!< 1D SPHERICAL
  CARTESIAN,    //!< 1D (slab) or 2D (xy) cartesian geometry
  END_GEOMETRY  //!< Sentinel value
};

} // end namespace rtt_mesh_element

#endif // mesh_element_Geometry_hh

//---------------------------------------------------------------------------//
// end of mesh_elementGeometry.hh
//---------------------------------------------------------------------------//
