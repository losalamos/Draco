//----------------------------------*-C++-*----------------------------------//
// AMR_Interface.cc
// Thomas M. Evans
// Thu Jul 16 09:25:43 1998
//---------------------------------------------------------------------------//
// @> AMR (Rage) interface implementation
//---------------------------------------------------------------------------//

#include "imc/AMR_Interface.hh"
#include "imc/AMR_Builder.hh"
#include "imc/Host_Manager.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

//===========================================================================//
// RAGE INTERFACE FUNCTION
//===========================================================================//

//---------------------------------------------------------------------------//
// F90 Rage interface functions
//---------------------------------------------------------------------------//
// IMC launcher in Rage

void rage_imc_(int *imc_ncycle, int *local_numtop, int *global_numtop,
	       double *imc_node_coord, int *imc_layout, int *num_b_cells,
	       int *b_proc, int *b_cell, double *imc_cv, double *imc_tev, 
	       double *imc_rev, double *imc_rho, double *imc_mua_n,
	       double *imc_mut_n, double *imc_dt, double *imc_time, 
	       double *imc_implicitness, int *imc_np_nom, int *imc_np_max, 
	       double *imc_dnpdt, int *imc_random_seed, 
	       int *imc_buffer_size, int *global_cells, double *imc_t4_slope,
	       double *return_e_dep, double *return_rad_den)
{
  // stl components
    using std::cout;
    using std::endl;
    using std::ofstream;
    using std::ostringstream;
    using std::string;

  // draco components
    using IMC::OS_Mesh;
    using IMC::AMR_Interface;
    using IMC::AMR_Builder;
    using IMC::Host_Manager;
    using C4::Wtime;
    using C4::node;
    using dsxx::SP;

  // timing info
    double begin;
    double end;

    if (node() == 0)
    {
      // welcome to IMC
	cout << "*********************************" << endl;
	cout << ">>> MILAGRO, 'a true miracle' <<<" << endl;
	cout << ">>> version 1.0               <<<" << endl;
	cout << ">>> Evans and Urbatsch        <<<" << endl;
	cout << "*********************************" << endl;

      // begining time
	begin = Wtime();
    }

  // barrier sync so everything starts at the same time
    C4::gsync();

  // make arguments struct
    AMR_Interface::Arguments arg(imc_node_coord, imc_layout, b_proc, b_cell,
				 global_cells, imc_cv, imc_rho, imc_mua_n,
				 imc_tev, imc_rev, imc_t4_slope,
				 *local_numtop,  *global_numtop,
				 *num_b_cells, *imc_implicitness, *imc_dt,  
				 *imc_time, *imc_dnpdt, *imc_np_nom,
				 *imc_np_max, 0.0, *imc_random_seed,
				 *imc_buffer_size, 0, *imc_ncycle,
				 return_e_dep, return_rad_den);

  // make a Rage manager and run IMC
    Host_Manager<OS_Mesh, AMR_Builder, AMR_Interface> rage_mgr(*imc_ncycle);
    rage_mgr.execute_IMC(arg);

  // ending message
    cout << "IMC is done for cycle " << *imc_ncycle << endl;
 
  // ending time
    C4::gsync();
    if (node() == 0)
    {
        end = Wtime();
        std::cout << std::endl << ">> Problem Timing" << std::endl;
        std::cout.precision(4);
        std::cout << " ** We ran for " << std::setw(15)
                  << std::setiosflags(std::ios::scientific)
                  << end-begin << " seconds" << std::endl;
    }
}

//===========================================================================//
// class AMR_Interface
//===========================================================================//

IMCSPACE

//---------------------------------------------------------------------------//
// Arguments constructor
//---------------------------------------------------------------------------//

AMR_Interface::Arguments::Arguments(const double *node_coord_, 
				    const int *layout_, 
				    const int *b_proc_, const int *b_cell_,
				    const int *global_cell_, 
				    const double *dedt_, const double *rho_, 
				    const double *opacity_abs_, 
				    const double *tev_, const double *rev_, 
				    const double *t4_slope_,
				    int num_cells_, int global_num_cells_,
				    int num_b_cells_, double implicitness_, 
				    double delta_t_, double elapsed_t_,
				    double dnpdt_, int npnom_, int npmax_, 
				    double rad_s_tend_, int seed_, 
				    int buffer_, int print_f_, int cycle_,
				    double * const e_dep_, 
				    double * const rad_den_)
    : node_coord(node_coord_), layout(layout_), b_proc(b_proc_),
      b_cell(b_cell_), global_cell(global_cell_), dedt(dedt_), rho(rho_),
      opacity_abs(opacity_abs_), tev(tev_), rev(rev_), t4_slope(t4_slope_),
      num_cells(num_cells_), global_num_cells(global_num_cells_),
      num_b_cells(num_b_cells_), implicitness(implicitness_), 
      delta_t(delta_t_), elapsed_t(elapsed_t_), dnpdt(dnpdt_), npnom(npnom_),
      npmax(npmax_), rad_s_tend(rad_s_tend_), seed(seed_), buffer(buffer_),
      print_f(print_f_), cycle(cycle_), e_dep(e_dep_), rad_den(rad_den_)
{
    Require (num_cells != 0);
    Require (num_cells <= global_num_cells);
}

//---------------------------------------------------------------------------//
// AMR_Interface Constructor
//---------------------------------------------------------------------------//

AMR_Interface::AMR_Interface(const Arguments &arg)
    : arguments(arg), evol_ext(arguments.num_cells, 0.0), ss_pos(0), 
      ss_temp(0)
{
    Require (arguments.num_cells != 0);
}

//---------------------------------------------------------------------------//
// static members holding global data needed between timesteps
//---------------------------------------------------------------------------//

SP<Particle_Buffer<Particle<OS_Mesh> >::Census> AMR_Interface::host_census;

//---------------------------------------------------------------------------//
// AMR_Interface accessor functions for Opacity Builder
//---------------------------------------------------------------------------//
// return a vector holding the densities

vector<double> AMR_Interface::get_density() const
{
    vector<double> rho(arguments.num_cells);
    for (int i = 0; i < arguments.num_cells; i++)
	rho[i] = arguments.rho[i];
    return rho;
}

//---------------------------------------------------------------------------//
// return an opacity vector (absorption opacities)

vector<double> AMR_Interface::get_kappa() const
{
    vector<double> kappa(arguments.num_cells);
    for (int i = 0; i < arguments.num_cells; i++)
	kappa[i] = arguments.opacity_abs[i];
    return kappa;
}

//---------------------------------------------------------------------------//
// return a Thompson scattering cross section (not getting any now)

vector<double> AMR_Interface::get_kappa_thomson() const
{
    vector<double> kappa(arguments.num_cells);
    for (int i = 0; i < arguments.num_cells; i++)
	kappa[i] = 0.0;
    return kappa;
}

//---------------------------------------------------------------------------//
// return a specific_heat (dE/dt) vector

vector<double> AMR_Interface::get_specific_heat() const
{
    vector<double> cv(arguments.num_cells);
    for (int i = 0; i < arguments.num_cells; i++)
	cv[i] = arguments.dedt[i];
    return cv;
}

//---------------------------------------------------------------------------//
// return temperature vectors

vector<double> AMR_Interface::get_temperature() const
{
    vector<double> T(arguments.num_cells);
    for (int i = 0; i < arguments.num_cells; i++)
	T[i] = arguments.tev[i];
    return T;
}

//---------------------------------------------------------------------------//
// AMR_Interface accessor functions for Source_Init
//---------------------------------------------------------------------------//
// return radiation temperature vector

vector<double> AMR_Interface::get_rad_temp() const
{
    vector<double> radT(arguments.num_cells);
    for (int i = 0; i < arguments.num_cells; i++)
	radT[i] = arguments.rev[i];
    return radT;
}

//---------------------------------------------------------------------------//
// return a global cell-list on the local cell map

vector<int> AMR_Interface::get_global_cells() const
{
    vector<int> global_cell(arguments.num_cells);
    for (int i = 0; i < arguments.num_cells; i++)
	global_cell[i] = arguments.global_cell[i];
    return global_cell;
}

//---------------------------------------------------------------------------//
// return the T^4 slope in each cell

vector<vector<double> > AMR_Interface::get_t4_slope() const
{
  // hardwired 3-D mesh
    vector<vector<double> > t4(3);
    t4[0].resize(arguments.num_cells);
    t4[1].resize(arguments.num_cells);
    t4[2].resize(arguments.num_cells);

  // loop through and assign the t4_slope
    int counter = 0;
    for (int cell = 0; cell < arguments.num_cells; cell++)
	for (int d = 0; d < t4.size(); d++)
	    t4[d][cell] = arguments.t4_slope[counter++];
    Check (counter == arguments.num_cells * t4.size());

  // return T^4 slope
    return t4;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of AMR_Interface.cc
//---------------------------------------------------------------------------//
