//----------------------------------*-C++-*----------------------------------//
// testP13T.hh
// Randy M. Roberts
// Fri Mar 20 12:04:07 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_testP13T_testP13T_hh__
#define __3T_testP13T_testP13T_hh__

#include "ds++/SP.hh"
#include "3T/P13T.hh"
#include "3T/testP13T/DiffusionSolverStub.hh"
#include "3T/testP13T/MeshTypeStub.hh"

// #include "3T/testP13T/MaterialPropertiesStub.hh"
#include "matprops/InterpedMaterialProps.hh"

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class testP13T - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class testP13T
{

    // NESTED CLASSES AND TYPEDEFS

    typedef MeshTypeStub MT;
    typedef InterpedMaterialProps MP;
    typedef DiffusionSolverStub<MT> DS;
    
    // typedef P13T<MT,MP,DS> P13T;
    
    // DATA

    SP< P13T<MT,MP,DS> > spP13T;
    
  public:

    // CREATORS
    
    testP13T();

    // MANIPULATORS
    
    // ACCESSORS

    void solve() const;

  private:
    
    // IMPLEMENTATION
};

END_NS_XTM  // namespace XTM

#endif                          // __3T_testP13T_testP13T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/testP13T.hh
//---------------------------------------------------------------------------//
