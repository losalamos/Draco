//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Surface_Sub_Tally.cc
 * \author Mike Buksas
 * \date   Mon Jun 23 15:33:15 2003
 * \brief  Implementation file for Surface_Sub_Tally
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Surface_Sub_Tally.hh"
#include "Azimuthal_Mesh.hh"
#include "Surface_Tracking_Interface.hh"

using rtt_dsxx::SP;
using std::vector;

namespace rtt_imc
{

//---------------------------------------------------------------------------//
/*! 
 * \brief Construct a surface sub-tally from the mesh and the global number
 * of surfaces
 * 
 * \param az_mesh The azimuthal mesh
 * \param surfaces_ The global number of surfaces
 *
 */
Surface_Sub_Tally::Surface_Sub_Tally(SP<Azimuthal_Mesh> az_mesh, int surfaces_)
    : azimuthal_mesh(az_mesh),
      num_tallies(2 * surfaces_),
      surfaces(surfaces_),
      weight_tally(2 * surfaces_),
      count_tally(2 * surfaces_),
      mesh_size( az_mesh->size() )
{

    Require(surfaces > 0)
    Require(azimuthal_mesh);

    // Initialize the size of each individual tally.
    for (int i = 0; i != num_tallies; ++i) 
    {
	weight_tally[i].resize(mesh_size, 0.0);
	count_tally[i].resize(mesh_size, 0);
    }

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Construct a surface sub-tally from the mesh and the global number
 * of surfaces
 * 
 * \param az_mesh The azimuthal mesh
 * \param surfaces_ The global number of surfaces
 *
 */
Surface_Sub_Tally::Surface_Sub_Tally(SP<Azimuthal_Mesh> az_mesh, 
				     const Surface_Tracking_Interface& interface)
    : azimuthal_mesh(az_mesh),
      mesh_size( az_mesh->size() )
{

    surfaces = interface.number_of_surfaces();
    num_tallies = 2 * surfaces;
    weight_tally.resize(num_tallies);
    count_tally.resize(num_tallies);

    Require(surfaces > 0)
    Require(azimuthal_mesh);

    // Initialize the size of each individual tally.
    for (int i = 0; i != num_tallies; ++i) 
    {
	weight_tally[i].resize(mesh_size, 0.0);
	count_tally[i].resize(mesh_size, 0);
    }

}

//---------------------------------------------------------------------------//
Surface_Sub_Tally::~Surface_Sub_Tally() { /* ... */ }


//---------------------------------------------------------------------------//
/*! 
 * \brief Add a surface crossing to the tally for the given surface
 *
 *
 * \param surface Which surface is crossed. 1-based.
 * \param direction Direction of the particle at crossing
 * \param is_outward True, if the particle is crossing outward
 * \param ew Energy-weight of the particle at the point of crossing.
 * \return void
 */
void Surface_Sub_Tally::add_to_tally(int surface,
				     const vector<double>& direction,
				     bool is_outward, double ew)
{

    Check ( ew > 0 );
    Check ( surface > 0 ); 
    Check ( surface <= surfaces );

    // Compute the index of the surface & direction in the tally.
    //    int surface_index = 2 * surface + static_cast<int>(is_outward);
    //    0-based surface_index must be < number of tallies
    int surface_index = get_surface_index(surface, is_outward);

    Check ( surface_index >= 0 ); 
    Check ( surface_index < num_tallies );

    // Get the bin from the mesh object
    int bin = azimuthal_mesh->find_bin(direction);
 
    Require ( bin >= 1 ); 
    Require ( bin <= mesh_size );

    // Compute the index of the bin in the tally.
    int bin_index = bin - 1;
    
    weight_tally[surface_index][bin_index] += ew;
    count_tally[surface_index][bin_index] += 1;

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the vector of outward weight data for a given surface.
 * 
 * \param surface The surface number. 1-based.
 * \return const vector<double> continaing the per-bin information.
 */
const std::vector<double>& Surface_Sub_Tally::get_outward_weight_tally(
    int surface) const
{
    Check (surface > 0);  Check(surface <= surfaces);

    int surface_index = 2 * (surface - 1) + 1;
    return weight_tally[surface_index];

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the vector of inward weight data for a given surface.
 * 
 * \param surface The surface number. 1-based.
 * \return const vector<double> containing the per-bin information.
 */
const std::vector<double>& Surface_Sub_Tally::get_inward_weight_tally(
    int surface) const
{
    Check (surface > 0);  Check(surface <= surfaces);

    int surface_index = 2 * (surface - 1) ;
    return weight_tally[surface_index];

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the vector of outward count data for a given surface.
 * 
 * \param surface The surface number. 1-based.
 * \return const vector<double> continaing the per-bin information.
 */
const std::vector<int>& Surface_Sub_Tally::get_outward_count_tally(int surface)
    const
{
    Check (surface > 0);  Check(surface <= surfaces);

    int surface_index = 2 * (surface - 1) + 1;
    return count_tally[surface_index];

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the vector of inward count data for a given surface.
 * 
 * \param surface The surface number. 1-based.
 * \return const vector<double> containing the per-bin information.
 */
const std::vector<int>& Surface_Sub_Tally::get_inward_count_tally(int surface)
    const
{
    Check (surface > 0);  Check(surface <= surfaces);

    int surface_index = 2 * (surface - 1) ;
    return count_tally[surface_index];

}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                 end of Surface_Sub_Tally.cc
//---------------------------------------------------------------------------//
