//----------------------------------*-C++-*----------------------------------//
// Diffusion_DB.cc
// Geoffrey M. Furnish
// Thu May 21 10:58:19 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "Diffusion_DB.hh"

#include "nml/Group.hh"
#include "nml/Items.hh"

namespace rtt_diffusion
{
 
void Diffusion_DB::setup_namelist( NML_Group& g )
{
#include ".nml_diffusion.cc"
}

} // end namespace rtt_diffusion

//---------------------------------------------------------------------------//
//                              end of Diffusion_DB.cc
//---------------------------------------------------------------------------//
