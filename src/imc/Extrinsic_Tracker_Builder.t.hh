//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Extrinsic_Tracker_Builder.t.hh
 * \author Mike Buksas
 * \date   Thu Jul 17 13:16:13 2003
 * \brief  Implementation file for Extrinsic_Tracker_Builder
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Extrinsic_Tracker_Builder_t_hh
#define rtt_imc_Extrinsic_Tracker_Builder_t_hh 

#include "mc/Sphere.hh"
#include "mc/Surface_Descriptor.hh"
#include "Extrinsic_Tracker_Builder.hh"

namespace rtt_imc
{

//---------------------------------------------------------------------------//
/*! 
 * \brief Construct an Extrinsic_Tracker_Builder using a RZWedge_Mesh and a
 * smart pointer to the Surface_Tracking_Interface
 * 
 * \param mesh_ The MT object
 * \param interface A Smart pointer to an implementation of
 * Surface_Tracking_Interface
 */
template<class MT>
Extrinsic_Tracker_Builder<MT>::Extrinsic_Tracker_Builder(
    const MT&        mesh_,
    const Mesh_Data& mesh_data_,
    SP_Interface     interface) 
    : mesh(mesh_),
      mesh_data(mesh_data_),
      number_of_cells(mesh_.num_cells()),
      global_surface_number(0),
      surfaces(),
      surface_indices(),
      surface_in_cell(mesh_.num_cells()),
      surface_areas()
{ 

    construction_implementation(*interface);
    
}
						     
//---------------------------------------------------------------------------//
/*! 
 * \brief Construct an Extrinsic_Tracker_Builder using an RZWedge_Mesh and a
 * Surface_Tracking_Interface
 * 
 * \param mesh_ The MT object
 * \param interface An implementation of Surface_Tracking_Interface
 */
template<class MT>
Extrinsic_Tracker_Builder<MT>::Extrinsic_Tracker_Builder(
    const MT&        mesh_,
    const Mesh_Data& mesh_data_,
    const Interface& interface)
    : mesh(mesh_),
      mesh_data(mesh_data_),
      number_of_cells(mesh_.num_cells()),
      global_surface_number(0),
      surfaces(),
      surface_indices(),
      surface_in_cell(mesh_.num_cells()),
      surface_areas()
{

    construction_implementation(interface);

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Builds and returns the surface tracker. 
 *
 * A tracker is built if there are any surfaces in the problem.  These
 * surfaces are defined globally.  A tracker on a given processor may not
 * have any surfaces; however, a tracker is still built.  The tracker may be
 * used to access the surface areas for \b all defined surfaces in the
 * problem.
 *
 * \return an instantiated tracker if there are any tally surfaces defined
 * globally in the problem.
 */
template<class MT>
typename Extrinsic_Tracker_Builder<MT>::SP_Tracker
Extrinsic_Tracker_Builder<MT>::build_tracker()
{
    // make a tracker
    SP_Tracker tracker;

    // build a tracker if surfaces are defined in the problem, during
    // contruction (in the constructor) the global_surface_number will be
    // equal to the number of global surfaces in the problem; if this is
    // non-zero then build the tracker
    if (global_surface_number) 
    {
	Check (global_surface_number == surface_areas.size());

	// build the tracker
	tracker = new Extrinsic_Surface_Tracker(global_surface_number,
						surfaces, 
						surface_indices,
						surface_areas,
						surface_in_cell);
	Ensure (tracker);
    }

    return tracker;
}

//---------------------------------------------------------------------------//
// Implementation:
//---------------------------------------------------------------------------//
/*!
 * \brief Implement the steps common to both constructors.
 */
template<class MT>
void Extrinsic_Tracker_Builder<MT>::construction_implementation(
    const Interface& interface)
{
    using rtt_mc::Surface_Descriptor;
    using std::vector;

    int given_surface_number = interface.number_of_surfaces();

    const vector<Surface_Descriptor>& surface_data(
	interface.get_surface_data());

    for (int surface = 0; surface < given_surface_number; ++surface)
    {
	global_surface_number++;
	process_surface( surface_data[surface] );
    }

    Ensure (global_surface_number == given_surface_number);
    Ensure (surface_areas.size() == global_surface_number);
    Ensure (surfaces.size() == surface_indices.size());
    Ensure (surfaces.size() <= given_surface_number);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a new surface.
 */
template<class MT>
void Extrinsic_Tracker_Builder<MT>::process_surface(
    const rtt_mc::Surface_Descriptor& descriptor)
{
    if (descriptor.type == rtt_mc::Surface_Descriptor::SPHERE)
	process_sphere(descriptor);
    else
	Insist(0, "Invalid surface descriptor encountered.");
}
		     
//---------------------------------------------------------------------------//

template<class MT>
void Extrinsic_Tracker_Builder<MT>::add_surface_to_list(SP_Surface surface)
{
    surfaces.push_back(surface);
    surface_indices.push_back(global_surface_number);
}

//---------------------------------------------------------------------------//

template<class MT>
bool Extrinsic_Tracker_Builder<MT>::check_point(
    const rtt_mc::Surface& surface,
    double                 x, 
    double                 z)
{
    Check(x >= 0);
    
    static std::vector<double> point(3, 0.0);

    point[0] = x;  point[1] = 0.0;  point[2] = z;

    return surface.is_inside(point);
}

//---------------------------------------------------------------------------//
// Sphere-specific functions:
//---------------------------------------------------------------------------//

template<class MT>
void Extrinsic_Tracker_Builder<MT>::process_sphere(
    const rtt_mc::Surface_Descriptor& d)
{
    
    double z(d.data[0]);
    double r(d.data[1]);

    Check ( r > 0);
    
    rtt_dsxx::SP<rtt_mc::Sphere> sphere ( new rtt_mc::Sphere(z, r) );

    bool on_mesh = check_intersections(*sphere);

    if (on_mesh) 
    {
	// if the sphere intersects a cell on this mesh then add it
	add_surface_to_list(sphere);
    }
    
    // build surface areas for all surfaces
    build_surface_areas(*sphere);
}

//---------------------------------------------------------------------------//

template<class MT>
bool Extrinsic_Tracker_Builder<MT>::check_intersections(
    const rtt_mc::Sphere& sphere)
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
 * \brief Generalization of sphere_intersects_cell.
 *
 * This function will always throw an assertion as specializations are
 * required to implement this.
 */
template<class MT>
bool Extrinsic_Tracker_Builder<MT>::sphere_intersects_cell(
    const rtt_mc::Sphere& sphere, 
    int                   cell)
{
    throw rtt_dsxx::assertion(
	"No specialization provided for this MT in Extrinsic_Tracker_Builder.");

    return false;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Generalization of build_surface_areas.
 *
 * This function will always throw an assertion as specializations are
 * required to implement this.
 */
template<class MT>
void Extrinsic_Tracker_Builder<MT>::build_surface_areas(
    const rtt_mc::Sphere& sphere)
{
    throw rtt_dsxx::assertion(
	"No specialization provided for this MT in Extrinsic_Tracker_Builder.");
}

} // end namespace rtt_imc

#endif // rtt_imc_Extrinsic_Tracker_Builder_t_hh

//---------------------------------------------------------------------------//
//                   end of imc/Extrinsic_Tracker_Builder.t.hh
//---------------------------------------------------------------------------//
