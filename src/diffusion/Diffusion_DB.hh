//----------------------------------*-C++-*----------------------------------//
// Diffusion_DB.hh
// Geoffrey M. Furnish
// Thu May 21 10:58:19 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __diffusion_Diffusion_DB_hh__
#define __diffusion_Diffusion_DB_hh__

class NML_Group;

namespace rtt_diffusion
{

//===========================================================================//
// class Diffusion_DB - 

// 
//===========================================================================//

class Diffusion_DB {
  public:

#include ".nml_diffusion.hh"

    void setup_namelist( NML_Group& g );
};

} // end namespace rtt_diffusion

#endif                          // __diffusion_Diffusion_DB_hh__

//---------------------------------------------------------------------------//
//                              end of diffusion/Diffusion_DB.hh
//---------------------------------------------------------------------------//
