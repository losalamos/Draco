//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Extrinsic_Tracker_Builder.cc
 * \author Mike Buksas
 * \date   Thu Jul 17 13:16:13 2003
 * \brief  Implementation file for Extrinsic_Tracker_Builder
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Extrinsic_Tracker_Builder.hh"
#include "mc/Sphere.hh"
#include "mc/RZWedge_Mesh.hh"

using rtt_mc::RZWedge_Mesh;
using rtt_mc::Sphere;
using rtt_dsxx::SP;
using std::vector;

namespace rtt_imc
{

Extrinsic_Tracker_Builder::Extrinsic_Tracker_Builder(const RZWedge_Mesh& mesh_) :
    mesh(mesh_), 
    number_of_cells(mesh_.num_cells()),
    global_surface_number(0),
    local_surfaces(0),
    surfaces(),
    surface_indices(),
    surface_in_cell(mesh_.num_cells())
{
    
    Check ( number_of_cells > 0);
    Check ( surface_in_cell.size() > 0);

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Builds and returns the surface tracker. Stops accumulation of
 * surfaces.
 */
    
SP<Extrinsic_Surface_Tracker> Extrinsic_Tracker_Builder::build_tracker()
{

    if (!tracker);
    {

	tracker = 
	    new Extrinsic_Surface_Tracker(surfaces, surface_indices, 
					  surface_in_cell
		);

    }

    Ensure (tracker);

    return tracker;

}


//---------------------------------------------------------------------------//
/*! 
 * \brief Adds a sphere to the list of surfaces. 
 * 
 * \param z z-coordinate of the sphere
 * \param r raidus of the sphere
 */
void Extrinsic_Tracker_Builder::add_sphere(double z, double r)
{
    
    Check ( r > 0 ); 
    
    Insist(!tracker, "Attempt to add surfaces to existing surface tracker");

    SP<Sphere> sphere ( new Sphere(z, r) );

    global_surface_number++;

    bool on_mesh = check_intersections(*sphere);

    if (on_mesh) add_sphere_to_list(sphere);


}


//---------------------------------------------------------------------------//
// Implementation:
//---------------------------------------------------------------------------//

void Extrinsic_Tracker_Builder::add_sphere_to_list(SP<Sphere> sphere)
{
    surfaces.push_back(sphere);

    surface_indices.push_back(global_surface_number);

}


//---------------------------------------------------------------------------//
bool Extrinsic_Tracker_Builder::check_intersections(const Sphere& sphere)
{

    bool on_mesh(false);

    for (int cell = 1; cell <= number_of_cells; ++cell)
    {
	
	if (sphere_intersects_cell(sphere, cell) )
	{
	    surface_in_cell[cell-1] = true;
	    on_mesh = true;
	}
	
    }

    return on_mesh;

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Detects intersections of a given sphere and the indicated cell of
 * the mesh.
 * 
 * \param sphere The sphere object
 * \param cell The cell of the RZWedge_Mesh
 * \return true if the cell is intersected by the sphere
 */
bool Extrinsic_Tracker_Builder::sphere_intersects_cell(
    const rtt_mc::Sphere& sphere, int cell)
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
bool Extrinsic_Tracker_Builder::check_point(const Sphere& sphere,
					    double x, double z)
{
    Check(x >= 0);
    
    static vector<double> point(3, 0.0);

    point[0] = x;  point[1] = 0.0;  point[2] = z;

    return sphere.is_inside(point);

}


} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                 end of Extrinsic_Tracker_Builder.cc
//---------------------------------------------------------------------------//
