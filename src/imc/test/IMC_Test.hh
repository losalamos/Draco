//----------------------------------*-C++-*----------------------------------//
// IMC_Test.hh
// Thomas M. Evans
// Tue Apr 27 10:53:05 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Some services we will need to test the imc packages
//---------------------------------------------------------------------------//

#ifndef __imc_test_IMC_Test_hh__
#define __imc_test_IMC_Test_hh__

#include "../Interface.hh"
#include "mc/Layout.hh"
#include "mc/XYCoord_sys.hh"
#include "mc/OS_Mesh.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <vector>
#include <string>

namespace rtt_imc_test
{

//===========================================================================//
// FAILURE LIMIT
//===========================================================================//

inline bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    return false;
}

//===========================================================================//
// INTERFACE CLASS
//===========================================================================//
// make an interface for a 6 cell mesh

class IMC_Interface : public rtt_imc::Interface
{
    typedef std::vector<int> vi;
    typedef std::vector<double> vd;
    typedef std::vector<std::string> vs;
    typedef std::vector<std::vector<double> > vvd;
   
  private:
    // Mesh data
    std::string coord;
    vvd         fine_edge;
    vs          bnd;
    
    // data for the Opacity and Mat_State
    vd          density;
    vd          kappa;
    vd          kappa_thomson;
    vd          temperature;
    vd          specific_heat;
    double      implicitness;
    double      delta_t;
    std::string analytic_opacity;
    std::string analytic_sp_heat;

  public:
    // constructor
    inline IMC_Interface();

    // public copy functions for Mesh
    std::string get_coordinates() const {return coord;}
    vvd get_fine_edge() const {return fine_edge;} 
    vs get_boundaries() const {return bnd;}
    
    // public copy functions for Opacity<MT>
    vd get_density() const {return density;}
    vd get_kappa() const {return kappa;}
    vd get_kappa_thomson() const {return kappa_thomson;}
    vd get_specific_heat() const {return specific_heat;}
    vd get_temperature() const {return temperature;}
    inline std::string get_analytic_opacity() const;
    inline std::string get_analytic_sp_heat() const;

    // accessor function to get implicitness factor (Fleck's alpha)
    double get_implicit() const { return implicitness; }

    // public copy functions for Source_Init<MT>
    double get_delta_t() const { return delta_t; }

    // source interfaces
    double get_rad_s_tend() const { return double(); }
    vd get_rad_temp() const { return vd(); }
    int get_npmax() const { return int(); }
    int get_npnom() const { return int(); }
    double get_dnpdt() const { return double(); }
};

// constructor
IMC_Interface::IMC_Interface()
    :  rtt_imc::Interface(), coord("xy"), fine_edge(2), bnd(4),
       density(6), kappa(6), kappa_thomson(6), temperature(6),
       specific_heat(6), implicitness(1.0), delta_t(.001),
       analytic_opacity("straight"), analytic_sp_heat("straight")
{
    // make the Mesh stuff

    // calculate the fine edges
    fine_edge[0].resize(4);
    fine_edge[1].resize(3);

    fine_edge[0][0] = -1;
    fine_edge[1][0] = -1;

    for (int i = 1; i < fine_edge[0].size(); i++)
	fine_edge[0][i] = fine_edge[0][i-1] + 1;
    for (int i = 1; i < fine_edge[1].size(); i++)
	fine_edge[1][i] = fine_edge[1][i-1] + 2;

    // calculate the boundaries
    for (int i = 1; i < 4; i++)
	bnd[i] = "vacuum";
    bnd[0] = "reflect";
    
    // make the Opacity and Mat_State stuff

    for (int i = 0; i < 3; i++)
    {
	// density
	density[i]   = 1.0;
	density[i+3] = 2.0;

	// kappa (in cm^2/g)
	kappa[i]     = .1;
	kappa[i+3]   = .01;
	
	// kappa thomson
	kappa_thomson[i]   = .5;
	kappa_thomson[i+3] = 0.0;

	// specific heat
	specific_heat[i]   = .1;
	specific_heat[i+3] = .2;

	// temperature
	temperature[i]   = 10;
	temperature[i+3] = 20;
    }
}


std::string IMC_Interface::get_analytic_opacity() const 
{ 
    return analytic_opacity;
}

std::string IMC_Interface::get_analytic_sp_heat() const 
{ 
    return analytic_sp_heat; 
}

} // end namespace rtt_imc_test

#endif                          // __imc_test_IMC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of imc/test/IMC_Test.hh
//---------------------------------------------------------------------------//
