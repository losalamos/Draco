//----------------------------------*-C++-*----------------------------------//
// AMR_Interface.cc
// Thomas M. Evans
// Thu Jul 16 09:25:43 1998
//---------------------------------------------------------------------------//
// @> AMR (Rage) interface implementation
//---------------------------------------------------------------------------//

#include "imc/AMR_Interface.hh"
#include "imc/OS_Mesh.hh"
#include "imc/OS_Builder.hh"
#include "imc/OS_Interface.hh"
#include "imc/IMC_Manager.hh"
#include <iostream>

//---------------------------------------------------------------------------//
// F90 Rage interface functions
//---------------------------------------------------------------------------//
// IMC launcher in Rage

void rage_IMC_()
{
  // stl components
    using std::cout;
    using std::endl;

  // draco components
    using IMC::OS_Mesh;
    using IMC::OS_Builder;
    using IMC::OS_Interface;
    using IMC::IMC_Manager;

  // welcome to IMC
    cout << "*********************************" << endl;
    cout << ">>> MILAGRO, 'a true miracle' <<<" << endl;
    cout << ">>> version 1.0               <<<" << endl;
    cout << ">>> Evans and Urbatsch        <<<" << endl;
    cout << "*********************************" << endl;
}

//---------------------------------------------------------------------------//
//                              end of AMR_Interface.cc
//---------------------------------------------------------------------------//
