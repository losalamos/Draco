//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/test/tstEnsight_Translator.cc
 * \author Thomas M. Evans
 * \date   Mon Jan 24 11:12:59 2000
 * \brief  Ensight_Translator test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Viz_Test.hh"
#include "../Ensight_Translator.hh"
#include "../Release.hh"
#include "ds++/Assert.hh"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;
using rtt_viz::Ensight_Translator;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_viz_test::fail(__LINE__);

//---------------------------------------------------------------------------//

void ensight_dump()
{
    // build an Ensight_Translator
    Ensight_Translator translator(1000);

    // dimensions
    int ncells   = 9; 
    int nvert    = 64; 
    int ndim     = 3;
    int ndata    = 1;
    int nhexvert = 8;
    int nrgn     = 1;
    
    // do an Ensight Dump
    vector<vector<int> >    ipar(ncells, vector<int>(nhexvert));
    vector<vector<double> > ens_vrtx_data(nvert, vector<double>(ndata));
    vector<vector<double> > ens_cell_data(ncells, vector<double>(ndata));
    vector<vector<double> > pt_coor(nvert, vector<double>(ndim));
  
    vector<int>    iel_type(ncells, 
			    rtt_viz::eight_node_hexahedron);
    vector<int>    rgn_index(ncells, 1);
    vector<string> ens_vdata_names(ndata, "Temp");
    vector<string> ens_cdata_names(ndata, "Velocity");
    vector<string> rgn_name(nrgn, "Problem");
    vector<int>    rgn_data(nrgn, 1);
    vector<int>    iproc(ncells, 0);

    string prefix   = "testproblem";
    int icycle      = 1;
    double time     = .01;
    double dt       = .01;
    string gd_wpath = ".";

    translator.ensight_dump(prefix, icycle, time, dt, gd_wpath,
			    ipar, iel_type, rgn_index, pt_coor,
			    ens_vrtx_data, ens_cell_data, ens_vdata_names,
			    ens_cdata_names, rgn_data, rgn_name, iproc);
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_viz::release() << endl; 
	    return 0;
	}

    try
    {
	// tests
	ensight_dump();
    }
    catch(dsxx::assertion &ass)
    {
	cout << "Dumbass you screwed up on " << ass.what() << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "**********************************************" << endl;
    if (passed) 
    {
        cout << "**** Ensight_Translator Self Test: PASSED ****" << endl;
    }
    cout <<     "**********************************************" << endl;
    cout << endl;
    
    cout << "Done testing Ensight_Translator." << endl;
}

//---------------------------------------------------------------------------//
//                              end of tstEnsight_Translator.cc
//---------------------------------------------------------------------------//
