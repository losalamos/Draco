//----------------------------------*-C++-*----------------------------------//
// Global_Buffer.cc
// Thomas M. Evans
// Wed Jun 17 10:21:19 1998
//---------------------------------------------------------------------------//
// @> Global_Buffer class implementation file
//---------------------------------------------------------------------------//

#include "imc/Global_Buffer.hh"
#include <iomanip>

IMCSPACE

using std::setiosflags;
using std::ios;
using std::endl;
using std::setw;

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//

template<class MT>
Global_Buffer<MT>::Global_Buffer(const MT &mesh, const Mat_State<MT>
				 &material, const Source_Init<MT> &source)
    : temperature(mesh.num_cells()), Cv(mesh.num_cells()),
      evol_net(mesh.num_cells()) 
{
    Require (mesh.num_cells() == material.num_cells());

  // assign the data
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
    {
	temperature[cell-1] = material.get_T(cell);
	Cv[cell-1]          = material.get_Cv(cell);
	evol_net[cell-1]    = source.get_evol_net(cell);
    }
}

//---------------------------------------------------------------------------//
// functions to update the global-mesh data
//---------------------------------------------------------------------------//
// update the material temperatures

template<class MT>
void Global_Buffer<MT>::update_T(const vector<double> &tally)
{
    Require (tally.size() == temperature.size());

  // update cell temperatures
    for (int i = 0; i < temperature.size(); i++)
    {
      // calculate net energy deposition from radiation to electrons
	double delta_E = tally[i] - evol_net[i];

      // calculate new electron temperature in cell
	temperature[i] += delta_E / Cv[i];
    }
}

//---------------------------------------------------------------------------//
// update a mat_state

template<class MT>
void Global_Buffer<MT>::update_Mat(Mat_State<MT> &mat) const
{
    Require (num_cells() == mat.num_cells());

  // update the temperatures in the Mat_State
    for (int cell = 1; cell <= num_cells(); cell++)
	mat.get_T(cell) = temperature[cell-1];
}

//---------------------------------------------------------------------------//
// print diagnostics
//---------------------------------------------------------------------------//
// output the Global Buffer

template<class MT>
void Global_Buffer<MT>::print(ostream &out) const
{
    out << setw(10) << setiosflags(ios::right) << "Cell"
	<< setw(20) << setiosflags(ios::right) << "Temperature"
	<< endl;
    out << "==============================" << endl;
    out.precision(4);
    for (int i = 0; i < num_cells(); i++)
    {	
	out << setw(10) << i+1 << setw(20) << setiosflags(ios::scientific)
	    << temperature[i] << endl;
    }
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Global_Buffer.cc
//---------------------------------------------------------------------------//
