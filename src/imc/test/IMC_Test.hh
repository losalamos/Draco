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
  private:
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
    : density(6), kappa(6), kappa_thomson(6), temperature(6),
      specific_heat(6), implicitness(1.0), delta_t(.001),
      analytic_opacity("straight"), analytic_sp_heat("straight")
{
    // make the stuff
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

//===========================================================================//
// 2D Mesh Objects
//===========================================================================//

//---------------------------------------------------------------------------//
// 2D Layout of a 3x3 mesh

inline rtt_mc::Layout build_2DLayout()
{
    using rtt_mc::Layout;

    // size of our 3x3 layout
    Layout layout = 6;
    
    // build the layout for our six cell problem
    for (int i = 1; i <= layout.num_cells(); i++)
    {
	// dimension cells
	layout.set_size(i, 4);

	// assign
	layout(i, 1) = i - 1;
	layout(i, 2) = i + 1;
	layout(i, 3) = i - 3;
	layout(i, 4) = i + 3;
    }

    // do layout boundaries
    for (int i = 1; i <= 3; i++)
	layout(i, 3) = 0;
    for (int i = 4; i <= 6; i++)
	layout(i, 4) = 0;
    layout(1, 1) = 1;
    layout(4, 1) = 4;
    layout(3, 2) = 0;
    layout(6, 2) = 0;

    return layout;
}    

//---------------------------------------------------------------------------//
// build 2D mesh (3x3)

inline dsxx::SP<rtt_mc::OS_Mesh> build_2DMesh()
{
    using rtt_mc::OS_Mesh;
    using rtt_mc::XYCoord_sys;
    using rtt_mc::Layout;
    using dsxx::SP;
    using namespace std;

    // build the Coord_sys
    SP<XYCoord_sys> coord;
    coord = new XYCoord_sys();

    // build the Layout
    Layout layout = build_2DLayout();

    // build the vertices
    vector<vector<double> > vert(2);
    for (int i = 0; i < 2; i++)
	vert[i].resize(12);
    vert[0][0] = -1;
    vert[1][0] = -1;
    for (int i = 1; i < 4; i ++)
    {
	vert[0][i] = vert[0][i-1] + 1;
	vert[1][i] = -1;
    }
    for (int i = 4; i < 8; i++)
    {
	vert[0][i] = vert[0][i-4];
	vert[1][i] = 1;
    }
    for (int i = 8; i < 12; i++)
    {
	vert[0][i] = vert[0][i-4];
	vert[1][i] = 3;
    }

    // build the cell pair
    vector<vector<int> > cp(6);
    for (int i = 0; i < 6; i++)
	cp[i].resize(4);

    for (int i = 1; i <= 3; i++)
	for (int j = 1; j <= 2; j++)
	{
	    int cell       = 1 + (i-1) + 3*(j-1);
	    int ref_vertex = 1 + (i-1) + 4*(j-1);
	    
	    cp[cell-1][0] = ref_vertex;
	    cp[cell-1][1] = ref_vertex + 1;
	    cp[cell-1][2] = ref_vertex + 4;
	    cp[cell-1][3] = ref_vertex + 5;
	}

    // now make the Mesh
    SP<OS_Mesh> mesh(new OS_Mesh(coord, layout, vert, cp));
    return mesh;
}

} // end namespace rtt_imc_test

#endif                          // __imc_test_IMC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of imc/test/IMC_Test.hh
//---------------------------------------------------------------------------//
