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
	       double *imc_dnpdt, int *imc_random_seed, int *imc_buffer_size)
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
				 imc_cv, imc_rho, imc_mua_n, imc_tev,
				 imc_rev, *local_numtop, *num_b_cells,
				 *imc_implicitness, *imc_dt, *imc_time,
				 *imc_dnpdt, *imc_np_nom, *imc_np_max, 0.0,
				 *imc_random_seed, *imc_buffer_size, 0,
				 *imc_ncycle, *global_numtop);

  // make a Rage manager and run IMC
    Host_Manager<OS_Mesh, AMR_Builder, AMR_Interface> rage_mgr(*imc_ncycle);
    rage_mgr.execute_IMC(arg);
 
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

AMR_Interface::Arguments::Arguments(const double *nc, const int *l, 
				    const int *bp, const int *bc, 
				    const double *de, const double *r,
				    const double *op_abs, const double *t,
				    const double *re, int numc,
				    int numbc, double imp, double dt,
				    double et, double dndt, int nnom, 
				    int nmax, double rst, int s, int b, 
				    int pf, int cycle_, int gnumc)
    : node_coord(nc), layout(l), b_proc(bp), b_cell(bc), dedt(de), rho(r),
      opacity_abs(op_abs), tev(t), rev(re), num_cells(numc),
      num_b_cells(numbc), implicitness(imp), delta_t(dt), elapsed_t(et),
      dnpdt(dndt), npnom(nnom), npmax(nmax), rad_s_tend(rst), seed(s), 
      buffer(b), print_f(pf), cycle(cycle_), global_num_cells(gnumc) 
{
    Require (num_cells != 0);
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

CSPACE

//---------------------------------------------------------------------------//
//                              end of AMR_Interface.cc
//---------------------------------------------------------------------------//
