//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Mat_State.hh
 * \author Thomas M. Evans
 * \date   Mon Mar  9 16:06:28 1998
 * \brief  Mat_State class header file
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Mat_State_hh__
#define __imc_Mat_State_hh__

#include "ds++/Assert.hh"
#include <iostream>

//===========================================================================//
/*!
  
 * \class Mat_State

 * \brief Material state definition for Fleck and Cumming's IMC.

 * The Mat_State class holds all of the required material state data required
 * to perform a transport calculation over one timestep.  Like the
 * rtt_imc::Opacity class, it is designed to be rebuilt each timestep as the
 * material properties change.  The data provided is:

 * \arg density in g/cc
 * \arg temperature in keV
 * \arg differential internal energy--heat capacity--(dE/dT) in Jerks/keV
 * \arg specific heat capacity in Jerks/g/keV
  
 */
// revision history:
// -----------------
// original
// 20-MAR-2001 : cleaned up source code files with documentation
// 
// 
//===========================================================================//

namespace rtt_imc 
{

template<class MT>
class Mat_State
{
  public:
    // Useful typedefs
    typedef std::string              std_string;
    typedef typename MT::CCSF_double ccsf_double;

  private:
    // Density in g/cc.
    ccsf_double density;
    
    // Material temperature in keV.
    ccsf_double temperature;

    // Specific heat capacities in Jerks/g/keV.
    ccsf_double spec_heat;

  public:
    // Constructor.
    Mat_State(const ccsf_double &, const ccsf_double &, const ccsf_double &); 
    
    // >>> ACCESSORS
    
    //! Get the density in a cell in g/cc.
    double get_rho(int cell) const { return density(cell); }
    
    //! Get the temperature in a cell in keV.
    double get_T(int cell) const { return temperature(cell); }

    //! Get the differential internal energy in a cell in Jerks/keV.
    inline double get_dedt(int cell) const;

    //! Get the specific heat in Jerks/g/keV.
    double get_spec_heat(int cell) const { return spec_heat(cell); }

    //! Get the number of cells in the mesh
    int num_cells() const { return density.get_Mesh().num_cells(); }

    // >>> DIAGNOSTICS

    //! Print out the material state.
    void print(std::ostream &, int) const;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
// Return the dE/dT in a cell

template<class MT>
double Mat_State<MT>::get_dedt(int cell) const
{
    Require (cell > 0 && cell <= num_cells());

    // calculate dEdT from the specific heat
    double dEdT = spec_heat(cell) * spec_heat.get_Mesh().volume(cell) * 
	density(cell);

    Check (dEdT >= 0.0);
    return dEdT;
}

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT>
std::ostream& operator<<(std::ostream &output, const Mat_State<MT> &object)
{
  // print out opacities for all cells
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;
    
    output << setw(8)  << setiosflags(ios::right) << "Cell" 
	   << setw(15) << "Density" << setw(15) << "Temp" 
	   << setw(15) << "dEdT" << setw(15) << "Cv" << endl;
    output << "--------------------------------------"
	   << "------------------------------" << endl;

    for (int i = 1; i <= object.num_cells(); i++)
	object.print(output, i);
    return output;
}

} // end namespace rtt_imc

#endif                          // __imc_Mat_State_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Mat_State.hh
//---------------------------------------------------------------------------//
