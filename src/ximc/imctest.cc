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
#include "imc/IMC_Manager.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

int main(int argc, char *argv[])
{
  // IMC namespace declarations
    using IMC::OS_Mesh;
    using IMC::OS_Builder;
    using IMC::OS_Interface;
    using IMC::IMC_Manager;

  // std namespace declarations
    using std::string;
    using std::vector;

  // init C4 stuff
    C4::Init(argc, argv);

  // time indicators
    double begin;
    double end;

  // starting time
    if (C4::node() == 0)
	begin = C4::Wtime();

  // welcome
    if (!C4::node())
    {
	std::cout << "**************************************" << std::endl;
	std::cout << ">>> Welcome to IMCTEST version 1.0 <<<" << std::endl;
	std::cout << ">>> T.M. Evans and T.J. Urbatsch   <<<" << std::endl;
	std::cout << ">>> Group XTM, LANL                <<<" << std::endl;
	std::cout << "**************************************" << std::endl;
	std::cout << std::endl;
    }

  // determine command line arguments
    vector<string> arguments(argc);
    bool verbose = false;
    int input;
    for (int i = 0; i < argc; i++)
	arguments[i] = argv[i];
    for (int i = 0; i < argc; i++)
    {
      // get title
	if (arguments[i] == "-i")
	    input = i+1;
      
      // get verbosity string
	if (arguments[i] == "-v")
	    verbose = true;
    }
    arguments.resize(0);
    
  // make a manager
    IMC_Manager<OS_Mesh, OS_Builder, OS_Interface> manager(verbose);

  // execute IMC
    manager.execute_IMC(argv[input]);

  // ending time
    C4::gsync();
    if (C4::node() == 0)
    {
	end = C4::Wtime();
	std::cout << std::endl << ">> Problem Timing" << std::endl;
	std::cout.precision(4);
	std::cout << " ** We ran for " << std::setw(15) 
		  << std::setiosflags(std::ios::scientific)
		  << end-begin << " seconds" << std::endl; 
    }

  // c4 end
    C4::Finalize();

  // everything ok
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of imctest.cc
//---------------------------------------------------------------------------//
