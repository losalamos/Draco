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

#include "ds++/Assert.hh"
#include <iostream>
#include <string>

namespace rtt_imc 
{

template<class MT>
class Mat_State
{

private:
  // data which defines the material state
    typename MT::CCSF_double density;
    typename MT::CCSF_double temperature;
    typename MT::CCSF_double dedt;
    typename MT::CCSF_double spec_heat;
    std::string analytic_sp_heat;

public:
  // inline constructors
    Mat_State(const typename MT::CCSF_double &, const typename
	      MT::CCSF_double &, const typename MT::CCSF_double &,
	      const typename MT::CCSF_double &, const std::string); 
    
  // public member functions

  // return values of material state data
    double& get_rho(int cell) { return density(cell); }
    double get_rho(int cell) const { return density(cell); }
    double& get_T(int cell) { return temperature(cell); }
    double get_T(int cell) const { return temperature(cell); }
    double& get_dedt(int cell) { return dedt(cell); }
    double get_dedt(int cell) const { return dedt(cell); }
    double get_spec_heat(int cell) const { return spec_heat(cell); }
    std::string get_analytic_sp_heat() const { return analytic_sp_heat; }

  // get the number of cells in the mesh
    inline int num_cells() const;

  // diagnostic functions (for printing)
    void print(std::ostream &, int) const;
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT>
std::ostream& operator<<(std::ostream &, const Mat_State<MT> &);

//---------------------------------------------------------------------------//
// inline functions for Mat_State
//---------------------------------------------------------------------------//
// constructors

template<class MT>
Mat_State<MT>::Mat_State(const typename MT::CCSF_double &density_, 
			 const typename MT::CCSF_double &temp_,
			 const typename MT::CCSF_double &dedt_,
			 const typename MT::CCSF_double &spec_heat_,
			 const std::string analytic_sp_heat_)
    : density(density_), temperature(temp_), dedt(dedt_), 
      spec_heat(spec_heat_), analytic_sp_heat(analytic_sp_heat_)
{
    Ensure (density.get_Mesh() == temperature.get_Mesh());
    Ensure (density.get_Mesh() == dedt.get_Mesh());
    Ensure (density.get_Mesh() == spec_heat.get_Mesh());
}

//---------------------------------------------------------------------------//
// return the num_cells

template<class MT>
inline int Mat_State<MT>::num_cells() const
{
  // return the number of cells
    Check (density.get_Mesh().num_cells() ==
	   temperature.get_Mesh().num_cells());
    return density.get_Mesh().num_cells();
}

} // end namespace rtt_imc

#endif                          // __imc_Mat_State_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Mat_State.hh
//---------------------------------------------------------------------------//
