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

class IMC_Interface
{
    typedef std::vector<std::vector<double> > vec_double;
    typedef std::vector<std::string> vec_string;
   
  private:
    // Mesh data
    std::string coord;
    vec_double fine_edge;
    vec_string bnd;
    
    // data for the Opacity and Mat_State
    std::vector<double> density;
    std::vector<double> kappa;
    std::vector<double> kappa_thomson;
    std::vector<double> temperature;
    std::vector<double> specific_heat;
    double              implicitness;
    double              delta_t;
    std::string         analytic_opacity;
    std::string         analytic_sp_heat;

  public:
    // constructor
    inline IMC_Interface();

    // public copy functions for Mesh
    std::string get_coordinates() const {return coord;}
    vec_double get_fine_edge() const {return fine_edge;} 
    vec_string get_boundaries() const {return bnd;}
    
    // public copy functions for Opacity<MT>
    std::vector<double> get_density() const {return density;}
    std::vector<double> get_kappa() const {return kappa;}
    std::vector<double> get_kappa_thomson() const {return kappa_thomson;}
    std::vector<double> get_specific_heat() const {return specific_heat;}
    std::vector<double> get_temperature() const {return temperature;}
    std::string get_analytic_opacity() const { return analytic_opacity; }
    std::string get_analytic_sp_heat() const { return analytic_sp_heat; }

    // accessor function to get implicitness factor (Fleck's alpha)
    double get_implicit() const { return implicitness; }

    // public copy functions for Source_Init<MT>
    double get_delta_t() const { return delta_t; }
};

// constructor
IMC_Interface::IMC_Interface()
    :  coord("xy"), fine_edge(2), bnd(4),
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

} // end namespace rtt_imc_test

#endif                          // __imc_test_IMC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of imc/test/IMC_Test.hh
//---------------------------------------------------------------------------//
