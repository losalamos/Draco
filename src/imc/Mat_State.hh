//----------------------------------*-C++-*----------------------------------//
// Mat_State.hh
// Thomas M. Evans
// Mon Mar  9 16:06:28 1998
//---------------------------------------------------------------------------//
// @> Mat_State class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Mat_State_hh__
#define __imctest_Mat_State_hh__

//===========================================================================//
// class Mat_State - 
//
// Purpose : defines the material state of the medium, density and temperature
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "Names.hh"
#include "OS_Mesh.hh"
#include "SP.hh"
#include <cassert>
#include <iostream>

IMCSPACE

using std::ostream;

template<class MT>
class Mat_State
{
private:
  // data which defines the material state
    typename MT::CCSF density;
    typename MT::CCSF temperature;

public:
  // explicit constructor
    explicit Mat_State(const typename MT::CCSF &density_, const typename
		       MT::CCSF &temperature_)
	: density(density_), temperature(temperature_)
    {}

  // public member functions
    double& Density(int cell) { return density(cell); }
    double Density(int cell) const { return density(cell); }
    double& Temperature(int cell) { return temperature(cell); }
    double Temperature(int cell) const { return temperature(cell); }
    int Num_cells() const
    {
	assert (density.Mesh().Num_cells() ==
		temperature.Mesh().Num_cells());
	return density.Mesh().Num_cells();
    }

  // diagnostic print
    void Print(int) const;
};

// overloaded operators
template<class MT>
ostream& operator<<(ostream &, const Mat_State<MT> &);

CSPACE

#endif                          // __imctest_Mat_State_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Mat_State.hh
//---------------------------------------------------------------------------//
