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

#include "ds++/Soft_Equivalence.hh"
#include "mc/Sphere.hh"
#include "mc/RZWedge_Mesh.hh"
#include "mc/Global_Mesh_Data.hh"
#include "Global.hh"
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

//---------------------------------------------------------------------------//
// MESH SPECIALIZATIONS ON build_surface_areas()
//---------------------------------------------------------------------------//
/*! 
 * \brief Determines the surface areas of the spheres that are subtended
 *  by the mesh; specialized for RZWedge_Mesh.
 * 
 * \param sphere The sphere object
 */
template<>
void Extrinsic_Tracker_Builder<rtt_mc::RZWedge_Mesh>::build_surface_areas(
    const rtt_mc::Sphere& sphere)
{
    using std::vector;
    using rtt_dsxx::soft_equiv;

    // sphere radius and center
    double r_s = sphere.get_radius();
    double z_s = sphere.get_center();

    // the extents are [0->low x][1->high x][2->low z][3->high z]
    vector<double> extents = mesh_data.get_spatial_extents();
    Check (extents.size() == 4);

    // make sure that the sphere does not extend out the high x side of the
    // mesh
    Insist (z_s + r_s <= extents[1],
	    "Tally sphere outside of mesh on high x side");

    // make sure that the sphere is (a) entirely contained between (-z,+z) or
    // (b) is a half-sphere, origin = -z or +z
    bool   half_sphere         = false;
    double surface_area_factor = 1.0;
    if (soft_equiv(z_s, extents[2]) || soft_equiv(z_s, extents[3]))
    {
	surface_area_factor = 0.5;
	half_sphere         = true;
    }
    else
    {
	Insist (z_s - r_s >= extents[2], 
		"Tally sphere outside of mesh on low z side");
	Insist (z_s + r_s <= extents[3], 
		"Tally sphere outside of mesh on high z side");	
    }

    // calculate the surface area (the surface_area_factor is set to 1/2 if
    // we have a half sphere)
    double surface_area = surface_area_factor * 
	2.0 * mesh.get_theta_radians() * r_s * r_s;

    // add to the surface area list
    surface_areas.push_back(surface_area);
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                 end of Extrinsic_Tracker_Builder.cc
//---------------------------------------------------------------------------//
