//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Surface_Tally.cc
 * \author Mike Buksas
 * \date   Mon Jun 23 15:33:15 2003
 * \brief  Contains tally information about surface crossings
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Surface_Tally.hh"
#include "Azimuthal_Mesh.hh"

using rtt_dsxx::SP;
using std::vector;

namespace rtt_imc
{

Surface_Tally::Surface_Tally(SP<Azimuthal_Mesh> az_mesh)
    : the_mesh(az_mesh), 
      mesh_size(az_mesh->size()),
      inward(az_mesh->size()),
      outward(az_mesh->size())
{

    Require(the_mesh);

}

void Surface_Tally::add_to_tally(const vector<double>& direction,
				 bool is_outward, double ew)
{

    Check ( ew > 0 );

    int bin = the_mesh->find_bin(direction);
    
    Require ( bin >= 1); Require (bin <= mesh_size);

    if (is_outward) outward[--bin] += ew;
    else inward[--bin] += ew;

}

Surface_Tally::~Surface_Tally() { /* ... */ }


} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                 end of Surface_Tally.cc
//---------------------------------------------------------------------------//
