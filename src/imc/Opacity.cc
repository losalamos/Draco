//----------------------------------*-C++-*----------------------------------//
// Opacity.cc
// Thomas M. Evans
// Fri Feb  6 13:52:29 1998
//---------------------------------------------------------------------------//
// @> Opacity class implementation file
//---------------------------------------------------------------------------//

#include "Opacity.hh"
#include "OS_Mesh.hh"
#include "Layout.hh"

IMCSPACE

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
template<class MT>
void Opacity<MT>::quickset(double sig_t, double sig_a, double sig_s)
{
    int num_cells = sigma.Mesh().Layout().getNum_cell();
  // fill CCSF array with data
    for (int i = 1; i <= num_cells; i++)
    {
	sigma(i)   = sig_t;
	sigma_a(i) = sig_a;
	sigma_s(i) = sig_s;
    }
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Opacity.cc
//---------------------------------------------------------------------------//
