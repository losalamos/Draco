//----------------------------------*-C++-*----------------------------------//
// testFullP13T_DB.hh
// Randy M. Roberts
// Fri Jun 12 11:35:48 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_testP13T_testFullP13T_DB_hh__
#define __3T_testP13T_testFullP13T_DB_hh__

#include "ds++/String.hh"

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

class NML_Group;

BEGIN_NS_XTM
    
enum unitsSelect { SI, AstroPhysical };

//===========================================================================//
// class testFullP13T_DB - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class testFullP13T_DB
{
  public:

#include ".nml_P13T.hh"

  public:
    void setup_namelist( NML_Group& g );
};

END_NS_XTM  // namespace XTM

#endif                          // __3T_testP13T_testFullP13T_DB_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/testFullP13T_DB.hh
//---------------------------------------------------------------------------//
