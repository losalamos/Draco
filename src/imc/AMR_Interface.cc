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
#include "c4/global.hh"
#include <iostream>
#include <iomanip>

//---------------------------------------------------------------------------//
// F90 Rage interface functions
//---------------------------------------------------------------------------//
// IMC launcher in Rage

void rage_imc_(int *i)
{
  // stl components
    using std::cout;
    using std::endl;

  // draco components
    using IMC::OS_Mesh;
    using IMC::OS_Builder;
    using IMC::OS_Interface;
    using IMC::IMC_Manager;
    using C4::Wtime;
    using C4::node;

  // timing info
    double begin;
    double end;

    if (node() == 0)
    {
      // welcome to IMC
	cout << "*********************************" << endl;
	cout << ">>> MILAGRO, 'a true miracle' <<<" << endl;
	cout << ">>> version 1.0               <<<" << endl;
	cout << ">>> Evans and Urbatsch        <<<" << endl;
	cout << "*********************************" << endl;

      // begining time
	begin = Wtime();
    }

  // make a manager
    IMC_Manager<OS_Mesh, OS_Builder, OS_Interface> manager;

  // execute IMC
    manager.execute_IMC("marshak");

  // ending time
    C4::gsync();
    if (node() == 0)
    {
        end = Wtime();
        std::cout << std::endl << ">> Problem Timing" << std::endl;
        std::cout.precision(4);
        std::cout << " ** We ran for " << std::setw(15)
                  << std::setiosflags(std::ios::scientific)
                  << end-begin << " seconds" << std::endl;
    }

}

//---------------------------------------------------------------------------//
//                              end of AMR_Interface.cc
//---------------------------------------------------------------------------//
