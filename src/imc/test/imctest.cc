//----------------------------------*-C++-*----------------------------------//
// imctest.cc
// Thomas M. Evans
// Mon Jun 15 18:47:51 1998
//---------------------------------------------------------------------------//
// @> test our standalone IMC code, IMCTEST
//---------------------------------------------------------------------------//

#include "imc/Names.hh"
#include "imc/OS_Mesh.hh"
#include "imc/OS_Builder.hh"
#include "imc/OS_Interface.hh"
#include "imc/IMC_Man.hh"
#include "c4/global.hh"

using IMC::OS_Mesh;
using IMC::OS_Builder;
using IMC::OS_Interface;
using IMC::IMC_Man;

int main(int argc, char *argv[])
{
  // init C4 stuff
    C4::Init(argc, argv);

  // make a manager
    IMC_Man<OS_Mesh, OS_Builder, OS_Interface> manager(true);

  // initialize on the host
    manager.host_init(argv[1]);

  // initialize the IMC processors
    manager.IMC_init();

  // c4 end
    C4::Finalize();

  // everything ok
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of imctest.cc
//---------------------------------------------------------------------------//
