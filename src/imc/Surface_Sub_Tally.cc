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

using rtt_dsxx::SP;
using std::vector;

namespace rtt_imc
{

Surface_Sub_Tally::Surface_Sub_Tally(SP<Azimuthal_Mesh> az_mesh, int surfaces_)
    : azimuthal_mesh(az_mesh),
      mesh_size(az_mesh->size()),
      surfaces(surfaces_),
      tallies(2 * surfaces_),
      tally(2 * surfaces_)
{

    Require(surfaces > 0)
    Require(azimuthal_mesh);

    // Initialize the size of each individual tally.
    for (int i = 0; i != tallies; ++i) tally[i].resize(mesh_size, 0.0);

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Add a surface crossing to the tally for the given surface
 *
 *
 * \param surface Which surface is crossed. Zero-based.
 * \param direction Direction of the particle at crossing
 * \param is_outward True, if the particle is crossing outward
 * \param ew Energy-weight of the particle at the point of crossing.
 * \return void
 */
void Surface_Sub_Tally::add_to_tally(int surface, const vector<double>& direction,
				     bool is_outward, double ew)
{

    Check ( ew > 0 );
    Check ( surface >= 0 ); Check ( surface < surfaces );

    // Compute the index of the surface & direction in the tally.
    int surface_index = 2 * surface + static_cast<int>(is_outward);

    // Get the bin from the mesh object
    int bin = azimuthal_mesh->find_bin(direction);

    // Compute the index of the bin in the tally.
    int bin_index = bin - 1;
    
    Require ( bin > 1); Require (bin <= mesh_size);

    tally[surface_index][bin_index] += ew;

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the vector of outward tally data for a given surface.
 * 
 * \param surface The surface number. 1-based.
 * \return const vector<double> continaing the per-bin information.
 */
const std::vector<double>& Surface_Sub_Tally::get_outward_tally(int surface) const
{
    Check (surface > 0);  Check(surface <= surfaces);

    int surface_index = 2 * (surface - 1) + 1;
    return tally[surface_index];

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the vector of inward tally data for a given surface.
 * 
 * \param surface The surface number. 1-based.
 * \return const vector<double> containing the per-bin information.
 */
const std::vector<double>& Surface_Sub_Tally::get_inward_tally(int surface) const
{
    Check (surface > 0);  Check(surface <= surfaces);

    int surface_index = 2 * (surface - 1) ;
    return tally[surface_index];

}

Surface_Sub_Tally::~Surface_Sub_Tally() { /* ... */ }


} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                 end of Surface_Sub_Tally.cc
//---------------------------------------------------------------------------//
