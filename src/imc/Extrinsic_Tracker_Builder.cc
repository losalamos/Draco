//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Extrinsic_Tracker_Builder.cc
 * \author Mike Buksas
 * \date   Thu Jul 17 13:16:13 2003
 * \brief  Specializations on MT in Extrinsic_Tracker_Builder.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/Sphere.hh"
#include "mc/RZWedge_Mesh.hh"
#include "Extrinsic_Tracker_Builder.hh"

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// MESH SPECIALIZATIONS ON sphere_intersects_cells()
//---------------------------------------------------------------------------//
/*! 
 * \brief Detects intersections of a given sphere and the indicated cell of
 * the mesh; specialized for RZWedge_Mesh.
 * 
 * \param sphere The sphere object
 * \param cell The cell of the RZWedge_Mesh
 * \return true if the cell is intersected by the sphere
 */
template<>
bool Extrinsic_Tracker_Builder<rtt_mc::RZWedge_Mesh>::sphere_intersects_cell(
    const rtt_mc::Sphere& sphere, 
    int                   cell)
{

    double r_s = sphere.get_radius();
    double z_s = sphere.get_center();

    double x_min = mesh.get_low_x(cell);
    double x_max = mesh.get_high_x(cell);

    double z_min = mesh.get_low_z(cell);
    double z_max = mesh.get_high_z(cell); 
    
    Check ( x_min < x_max); Check (z_min < z_max);

    // Check if entire cell is left or right of sphere:
    if ( (z_max < z_s - r_s) || (z_min > z_s + r_s) )
    {
	return false;
    }
    else if ( (z_min < z_s) && (z_max > z_s) )
    { // The cell straddles the center
	
	if ( x_min > r_s ) return false;      // Cell is above the sphere.
	else if ( x_max > r_s) return true;   // Cell top is < sphere raidus.
	
	Require (x_max <= r_s);

	// Check for a "bulge" crossing, where x_max <= r_s but the corners
	// of the cell "stick out" of the sphere. If either of the upper
	// corners is outside the sphere, we have a bulge crossing.
	
	bool upper_left_inside  = check_point(sphere, x_max, z_min);
	bool upper_right_inside = check_point(sphere, x_max, z_max);

	return ((upper_left_inside == false) || (upper_right_inside == false)) ;
	
    }
    else if (z_max <= z_s) // The cell is to the left of the sphere.
    {
	// Check the status of the upper-left and lower-right corners. If
	// different, we have a crossing
	
	bool upper_left_inside  = check_point(sphere, x_max, z_min);
	bool lower_right_inside = check_point(sphere, x_min, z_max);
	
	return (upper_left_inside != lower_right_inside);
	
    }
    else if (z_min >= z_s) // To the right
    {
	// Check the status of the lower-left and upper-right corners. If
	// different, we have a crossing.
	
	bool lower_left_inside  = check_point(sphere, x_min, z_min);
	bool upper_right_inside = check_point(sphere, x_max, z_max);

	return (lower_left_inside != upper_right_inside);
	
    }
    else 
    {
	Insist(0, "Logic error in sphere_intersects_cell");
    }
    
    // default return condition
    return false;
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                 end of Extrinsic_Tracker_Builder.cc
//---------------------------------------------------------------------------//
