//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfInterpreter.cc
 * \author Allan Wollaber
 * \date   Fri Oct 12 15:39:39 2001
 * \brief  Basic reader to print info in IPCRESS files.
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "GandolfFile.hh"
#include "GandolfException.hh"
#include "GandolfGrayOpacity.hh"
#include "GandolfMultigroupOpacity.hh"
#include "cdi/OpacityCommon.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>

using rtt_dsxx::SP;
using rtt_cdi_gandolf::GandolfFile;
using rtt_cdi_gandolf::GandolfMultigroupOpacity;
using rtt_cdi_gandolf::GandolfGrayOpacity;
using rtt_cdi::GrayOpacity;
using rtt_cdi::MultigroupOpacity;
using std::cout;
using std::endl;
using std::string;
using std::ios;

//---------------------------------------------------------------------------//
/*!
 * \brief Basic reader to print info in IPCRESS files.
 *
 * Information about the GANDOLF IPCRESS Interface Library may be
 * found on X-Division's Physical Data Team web site:
 *
 * http://velvet.lanl.gov/PROJECTS/DATA/atomic/gandolf/intro.lasso
 *
 * We have a slightly modified copy of the gandolf libraries
 * (explicitly added compiler libraries) located at:
 *
 * /radtran/vendors/gandolf
 *
 * To use this package Draco must be compiled with Gandolf
 *
 * Modification of this executable to make it more useful is encouraged.
 */
//---------------------------------------------------------------------------//
void gandolf_file_read(const std::string& op_data_file) throw()
{
    // Gandolf data filename (IPCRESS format required)
	    
    cout << "Creating a Gandolf File object from " << op_data_file << endl;
	    
    SP<GandolfFile> spGF;
    try
    {
	spGF = new rtt_cdi_gandolf::GandolfFile( op_data_file );
    }
    catch ( const rtt_cdi_gandolf::gmatidsException& GandError )
    {
	std::cerr << "Error: Can't open file " << op_data_file << ". Aborting" << endl;
        throw;
	return;
    }

    // Test the new object to verify the constructor and accessors.
    size_t numMaterials = spGF->getNumMaterials();
    cout << "This opacity file has " << numMaterials << " materials:" << endl;

    // Print the Mat IDs.
    std::vector<int> matIDs = spGF->getMatIDs();
    for (size_t i=0; i<numMaterials; ++i)
        cout << "Material " << i+1 << " has ID number " << matIDs[i] << endl;

    size_t    i=0;
    cout << "=========================================" << endl;
    cout << "For Material " << i+1 << " with ID number " << matIDs[i] << endl;
    cout << "=========================================" << endl;
    SP<MultigroupOpacity> spMgOp;
    try
    {
       spMgOp = new GandolfMultigroupOpacity( spGF, matIDs[i], rtt_cdi::ROSSELAND, rtt_cdi::TOTAL);
    }
    catch ( rtt_cdi_gandolf::GandolfException& GandError)
    {
       std::cerr << "Failed to create gray GandolOpacity object for material " << 
                     matIDs[i] << endl << GandError.what();
    }

    // set precision
    cout.precision(7);
    cout.setf(ios::scientific);

    // Print the density grid
    std::vector<double> dens = spMgOp->getDensityGrid();
    cout   << "Density grid" << endl;
    for (size_t tIndex=0; tIndex < dens.size(); ++tIndex)
        cout << tIndex+1<< "\t" << dens[tIndex] << endl;

    // Print the temperature grid
    std::vector<double> temps = spMgOp->getTemperatureGrid();
    cout   << "Temperature grid" << endl;
    for (size_t tIndex=0; tIndex < temps.size(); ++tIndex)
        cout << tIndex+1 << "\t" << temps[tIndex] << endl;
    
    // Print the frequency grid
    std::vector<double> groups = spMgOp->getGroupBoundaries();
    cout   << "Frequency grid" << endl;
    for (size_t gIndex=0; gIndex < groups.size(); ++gIndex)
        cout << gIndex + 1<< "\t" << groups[gIndex] << endl;
 
     
    cout << endl;
 
   
    int keepGoing = 1;
    size_t matID(1);
    size_t denID(1);
    size_t tempID(1);
    size_t selID(1);
    while ( keepGoing )
    {
       cout << "Enter 0 to quit." << endl;
       cout << "Please select a material from 1 to " << numMaterials << endl;
       std::cin  >> keepGoing; if (keepGoing == 0) break;
       matID = keepGoing-1;  Insist( (  matID <numMaterials)  ,"Invalid material index"); 

       cout << "Please select a density from 1 to " << dens.size() << endl;
       std::cin  >> keepGoing; if (keepGoing == 0) break;
       denID = keepGoing-1;  Insist( ( denID <dens.size())  ,"Invalid density index"); 

       cout << "Please select a temperature from 1 to " << temps.size() << endl;
       std::cin  >> keepGoing; if (keepGoing == 0) break;
       tempID = keepGoing-1;  Insist( ( tempID < temps.size())  ,"Invalid temperature index"); 
 
       cout << "Choose your opacity type:" << endl;
       cout << "1: Rosseland Absorption (Gray), 2: Planck Absorption (Gray)" << endl;
       cout << "3: Rosseland Absorption (MG),   4: Planck Absorption (MG)"   << endl;
       std::cin  >> keepGoing; if (keepGoing == 0) break;
       selID = keepGoing-1;  Insist( ( selID < 4)  ,"Invalid choice."); 

       if (selID == 0)
       {
          SP<GrayOpacity> spGOp;
          spGOp = new GandolfGrayOpacity( spGF, matIDs[matID], rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
          cout << "The Gray Rosseland Absorption Opacity for " << endl;
          cout << "material " << matID << " Id(" << matIDs[matID] << ") at density "
                    << dens[denID] << ", temperature " << temps[tempID] << " is " 
                    << spGOp->getOpacity(temps[tempID], dens[denID]) << endl;
       }
       else if (selID == 1)
       {
          SP<GrayOpacity> spGOp;
          spGOp = new GandolfGrayOpacity( spGF, matIDs[matID], rtt_cdi::PLANCK, rtt_cdi::ABSORPTION);
          cout << "The Gray Planck Absorption Opacity for " << endl;
          cout << "material " << matID << " Id(" << matIDs[matID] << ") at density "
                    << dens[denID] << ", temperature " << temps[tempID] << " is " 
                    << spGOp->getOpacity(temps[tempID], dens[denID]) << endl;
          
       }
       else if (selID == 2)
       {
          SP<MultigroupOpacity> spMGOp;
          spMGOp = new GandolfMultigroupOpacity( spGF, matIDs[matID], rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
          cout << "The Multigroup Rosseland Absorption Opacity for " << endl;
          cout << "material " << matID << " Id(" << matIDs[matID] << ") at density "
                    << dens[denID] << ", temperature " << temps[tempID] << " is: " << endl ;
          std::vector<double> opData = spMGOp->getOpacity(temps[tempID], dens[denID]);
          cout << "Index \t Group Center \t\t Opacity" << endl;
          for (size_t g=0; g < opData.size(); ++g)
                  cout << g+1 << "\t " << 0.5*(groups[g]+groups[g+1]) << "   \t " << opData[g] << endl; 
       }
       else if (selID == 3)
       {
          SP<MultigroupOpacity> spMGOp;
          spMGOp = new GandolfMultigroupOpacity( spGF, matIDs[matID], rtt_cdi::PLANCK, rtt_cdi::ABSORPTION);
          cout << "The Multigroup Planck Absorption Opacity for " << endl;
          cout << "material " << matID << " Id(" << matIDs[matID] << ") at density "
                    << dens[denID] << ", temperature " << temps[tempID] << " is: " << endl ;
          std::vector<double> opData = spMGOp->getOpacity(temps[tempID], dens[denID]);
          cout << "Index \t Group Center  \t\t Opacity" << endl;
          for (size_t g=0; g < opData.size(); ++g)
                  cout << g+1 << "\t " << 0.5*(groups[g]+groups[g+1]) << "   \t " << opData[g] << endl; 
       }
    }
    
    cout << "Ending session." << endl; 
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    string filename;
    // Parse the arguments
    for (int arg = 1; arg < argc; arg++)
    {
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_dsxx::release() 
		 << endl;
	    return 0;
	}
	else if (string(argv[arg]) == "--help" || string(argv[arg]) == "-h")
	{
	    cout << argv[0] << ": version " << rtt_dsxx::release() << endl;
	    cout << "Usage: GandolfInterpreter <ipcress file>" << endl;
            cout << "Follow the prompts to print opacity data to the screen." << endl;
	    return 0;
	}
        else
            filename = string(argv[arg]);
    }

    try
    {
	// >>> UNIT TESTS
	gandolf_file_read(filename);
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While attempting to read an opacity file, " << ass.what()
	     << endl;
	return 1;
    }

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of GandolfInterpreter.cc
//---------------------------------------------------------------------------//
