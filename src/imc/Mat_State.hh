//----------------------------------*-C++-*----------------------------------//
// Mat_State.hh
// Thomas M. Evans
// Mon Mar  9 16:06:28 1998
//---------------------------------------------------------------------------//
// @> Mat_State class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Mat_State_hh__
#define __imc_Mat_State_hh__

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

#include "imc/Names.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <iostream>

IMCSPACE

using std::ostream;

template<class MT>
class Mat_State
{
public:
  // structure for holding the global temperatures between cycles
    struct Shell
    {
	vector<double> temperature;
	inline Shell(const Mat_State<MT> &);
    };

private:
  // data which defines the material state
    typename MT::CCSF_double density;
    typename MT::CCSF_double temperature;

public:
  // explicit constructor
    Mat_State(const typename MT::CCSF_double &density_, 
	      const typename MT::CCSF_double &temp_)
	: density(density_), temperature(temp_) {}

  // public member functions

  // return values of material state data
    double& get_rho(int cell) { return density(cell); }
    double get_rho(int cell) const { return density(cell); }
    double& get_T(int cell) { return temperature(cell); }
    double get_T(int cell) const { return temperature(cell); }

  // get the number of cells in the mesh
    inline int num_cells() const;

  // diagnostic functions (for printing)
    void print(ostream &, int) const;
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT>
ostream& operator<<(ostream &, const Mat_State<MT> &);

//---------------------------------------------------------------------------//
// inline functions for Mat_State<MT>::Shell
//---------------------------------------------------------------------------//

template<class MT>
inline Mat_State<MT>::Shell::Shell(const Mat_State<MT> &material)
    : temperature(material.num_cells())
{
    Require (temperature.size() == material.num_cells());

  // assign the temperatures to this vector array
    for (int cell = 1; cell <= temperature.size(); cell++)
	temperature[cell-1] = material.get_T(cell);
}

//---------------------------------------------------------------------------//
// inline functions for Mat_State
//---------------------------------------------------------------------------//

template<class MT>
inline int Mat_State<MT>::num_cells() const
{
  // return the number of cells
    Check (density.get_Mesh().num_cells() ==
	   temperature.get_Mesh().num_cells());
    return density.get_Mesh().num_cells();
}

CSPACE

#endif                          // __imc_Mat_State_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Mat_State.hh
//---------------------------------------------------------------------------//
