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
#include "Constants.hh"

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

 * \param beta_degrees_ "angle" of pyramid in degrees

*/

Pyramid_Mesh::Pyramid_Mesh(SP_Coord coord_,
			   AMR_Layout &layout_,
			   vf_double & cell_x_extents_,
			   double beta_degrees_)
    :coord(coord_),
     layout(layout_),
     cell_x_extents(cell_x_extents_),
     beta_degrees(beta_degrees_)
{
    // Check coordinate system class
    Require (coord);
    Require (coord->get_dim() == 3);
    Require (coord->get_Coord() == std:string("xyz"));

    // Make sure that the pyramid angle is postive and not obtuse
    Require ((beta_degrees>0.0) && (beta_degrees <= 90.0));

    // check that the cell-extents vector has num_cells elements
    // and weakly check that each cell-element has 2 extent-elements
    Check (cell_x_extents.size() == layout.num.cells());
    Check (cell_x_extents[0].size() == 2);
    Check (cell_x_extents[layout.num_cells()-1].size() ==2);

    // precalculate heavily used angle quatities
    calc_angle_data(beta_degrees);

    // calculate and assign the on-processor total volume
    //calc_total_volume();
}

 
//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate wedge angle data once for use throughout calculation
 *
 * \param beta_degrees_ "angle" of pyramid in degrees
 */
void Pyramid_Mesh::calc_angle_data(const double beta_degrees)
{
    Require ((beta_degrees >0.0) && (beta_degrees <= 90.0));
    
    beta_radians=beta_degrees*rtt_mc::global::pi/180.0;
    Check (beta_radians>0.0);
    Check (beta_radians <=rtt_mc::global::pi/2.0);

    // precalculate heavily used trif functions
    tan_beta=std::tan(beta_radians);
    sin_beta=std::sin(beta_radians);
    cos_beta=std::cos(beta_radians);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate and set the total (on-processor) volume of the mesh.
 *
 * \return total volume of on-processor Pyramid cells
 */
void Pyramid_Mesh::calc_total_volume()
{
    Require (num_cells()>0);

    // initialize private data member
    total_volume=0.0;
    
    // sum local cell volumes
    for(int cell=1; cell<=num_cells(); cell++)
	total_volume+=volume(cell);

    Ensure (total_volume >0.0);
}

} // end namespace rtt_mc

#endif

//---------------------------------------------------------------------------//
//                 end of Pyramid_Mesh.cc
//---------------------------------------------------------------------------//
