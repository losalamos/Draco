//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Pyramid_Mesh.cc
 * \author Jeffery Densmore
 * \date   Mon Oct  6 09:15:12 2003
 * \brief  Pyramid_Mesh implementation file.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Pyramid_Mesh_cc
#define rtt_mc_Pyramid_Mesch_cc

#include "Pyramid_Mesh.hh"
#include "dx++/Assert.hh"

namespace rtt_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*! 
 * \brief Pyramid_Mesh constructor.

 * The constructor requires complete arguments of the constituent data
 * necessary to build a Pyramid_Mesh.  It is expected that an appropriate
 * builder class will message general data through an interface to build the 
 * specific data structures needed by Pyramid_Mesh

 * \param coord_ coordinate system smart pointer

 * \param layout_ AMR_layout giving cell connectivity information

 * \param cell_x_extents_ the x coordinates of each cell in the following
 * form: [cell][lowx]; [cell][hi x];

 * \param alpha_degrees_ angle of the original spherical cone in degrees

 * \param beta_degrees_ "angle" of pyramid in degrees

*/

Pyramid_Mesh::Pyramid_Mesh(SP_Coord coord_,
			   AMR_Layout &layout_,
			   vf_double & cell_x_extents_,
			   double alpha_degrees_,
			   double beta_degrees_)
    :coord(coord_),
     layout(layout_),
     cell_x_extents(cell_x_extents_),
     alpha_degrees(alpha_degrees_),
     beta_degrees(beta_degrees_)
{
    // Check coordinate system class
    Require (coord);
    Require (coord->get_dim() == 3);
    Require (coord->get_Coord() == std:string("xyz"));

    // Make sure that both angles are postive and not obtuse
    Require ((alpha_degrees> 0.0) && (alpha_degrees <= 90.0));
    Require ((beta_degrees>0.0) && (beta_degrees <= 90.0));

    // check that the cell-extents vector has num_cells elements
    // and weakly check that each cell-element has 2 extent-elements
    Check (cell_x_extents.size() == layout.num.cells());
    Check (cell_x_extents[0].size() == 2);
    Check (cell_x_extents[layout.num_cells()-1].size() ==2);

    // precalculate heavily used angle quatities
    //calc_wedge_angle_data(alpha_degrees);

    // calculate and assign the on-processor total volume
    //calc_total_volume();
}

 




} // end namespace rtt_mc

#endif

//---------------------------------------------------------------------------//
//                 end of Pyramid_Mesh.cc
//---------------------------------------------------------------------------//
