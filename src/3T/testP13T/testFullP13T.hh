//----------------------------------*-C++-*----------------------------------//
// testFullP13T.hh
// Randy M. Roberts
// Sat May 30 15:09:58 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_testP13T_testFullP13T_hh__
#define __3T_testP13T_testFullP13T_hh__

#include "testFullP13T_DB.hh"
#include "ds++/SP.hh"
#include "3T/P13T.hh"
#include "3T/Diffusion_P1.hh"
#include "mesh/Mesh_XYZ.hh"
#include "matprops/InterpedMaterialProps.hh"

#include <string>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM

//===========================================================================//
// class testFullP13T - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class testFullP13T
{

    // NESTED CLASSES AND TYPEDEFS

  public:
    
    typedef Mesh_XYZ MT;
    typedef InterpedMaterialProps MP;
    typedef Diffusion_P1<MT> DS;

    // DATA

  private:
    
    SP< P13T<MT,MP,DS> > spP13T;

    testFullP13T_DB pdb;
    Diffusion_DB diffdb;
    pcg_DB pcg_db;
    
  public:

    // CREATORS
    
    testFullP13T(const std::string &infile);

    // MANIPULATORS
    
    // ACCESSORS

    void run() const;

  private:
    
    // IMPLEMENTATION
};

END_NS_XTM  // namespace XTM

#endif                          // __3T_testP13T_testFullP13T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/testFullP13T.hh
//---------------------------------------------------------------------------//
