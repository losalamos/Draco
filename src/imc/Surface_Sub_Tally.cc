//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Surface_Sub_Tally.cc
 * \author Mike Buksas
 * \date   Mon Jun 23 15:33:15 2003
 * \brief  Contains tally information about surface crossings
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
    : the_mesh(az_mesh),
      mesh_size(az_mesh->size()),
      surfaces(surfaces_),
      tallies(2 * surfaces_),
      tally(2 * surfaces_)
{

    Require(surfaces > 0)
    Require(the_mesh);

    for (int i = 0; i != tallies; ++i) tally[i].resize(mesh_size, 0.0);

}

void Surface_Sub_Tally::add_to_tally(int surface, const vector<double>& direction,
				 bool is_outward, double ew)
{

    Check ( ew > 0 );
    Check ( surface >= 0 ); Check ( surface < surfaces );

    int surface_index = 2 * surface + static_cast<int>(is_outward);
    int bin = the_mesh->find_bin(direction);
    int bin_index = bin - 1;
    
    Require ( bin > 1); Require (bin <= mesh_size);

    tally[surface_index][bin_index] += ew;

}

const vector<double>& Surface_Sub_Tally::get_outward_tally(int surface) const
{
    Check (surface > 0);  Check(surface <= surfaces);

    int surface_index = 2 * (surface - 1) + 1;
    return tally[surface_index];

}

const vector<double>& Surface_Sub_Tally::get_inward_tally(int surface) const
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
