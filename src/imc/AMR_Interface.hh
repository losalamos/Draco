//----------------------------------*-C++-*----------------------------------//
// AMR_Interface.hh
// Thomas M. Evans
// Thu Jul 16 09:25:43 1998
//---------------------------------------------------------------------------//
// @> Rage Interface functions and classes
//---------------------------------------------------------------------------//

#ifndef __imc_AMR_Interface_hh__
#define __imc_AMR_Interface_hh__

//===========================================================================//
// class AMR_Interface - 
//
// Purpose : Interface functions to Rage's AMR-Eulerian Code.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imc/Names.hh"
#include "imc/OS_Mesh.hh"
#include "imc/Particle.hh"
#include "imc/Particle_Buffer.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

//===========================================================================//
// F90 Functional Interface to Rage 
//===========================================================================//
// direct functional interface to F90 Rage code

extern "C"
{
    extern void milstone_(int *, int *, int *, double *, int *, int *, int *, 
			  int *, double *, double *, double *, double *, 
			  double *, double *, double *, double *, double *, 
			  int *, int *, double *, int *, int *, int *, 
			  double *, double *, double *);
}

//---------------------------------------------------------------------------//
// census object declarations

IMCSPACE

// stl components
using std::vector;
using std::string;

// draco components
using dsxx::SP;

//===========================================================================//
// class AMR_Interface
//===========================================================================//

class AMR_Interface
{
public:
  // structure for passing arguments to the interface, we do not reassign
  // memory here in copy constructors or assignment operators on purpose, our 
  // goal is to conserve some memory as we know that the initial Arguments
  // object made in the rage_imc function WILL NOT be destroyed during a
  // timestep
    struct Arguments
    {
      // data determining the problem layout

      // arrays that we receive
	const double *node_coord;
	const int *layout;
	const int *b_proc;
	const int *b_cell;
	const int *global_cell;
	const double *dedt;
	const double *rho;
	const double *opacity_abs;
	const double *tev;
	const double *rev;
	const double *t4_slope;
	
      // variables
	const int num_cells;
	const int global_num_cells;
	const int num_b_cells;	
	const double implicitness;
	const double delta_t;
	const double elapsed_t;
	const double dnpdt;
	const int npnom;
	const int npmax;
	const double rad_s_tend;
	const int seed;
	const int buffer;
	const int print_f;
	const int cycle;

      // arrays that we return data to
	double * const e_dep;
	double * const rad_den;

      // constructor
	Arguments(const double *, const int *, const int *, const int *,
		  const int *, const double *, const double *, 
		  const double *, const double *, const double *, 
		  const double *, int, int, int, double, double, double, 
		  double, int, int, double, int, int, int, int, 
		  double * const, double * const);
    };

private:
  // data from Rage
    Arguments arguments;

  // blank vector argument to be filled in at later date when the proper
  // features are added
    vector<double> evol_ext;
    vector<string> ss_pos;
    vector<double> ss_temp;

  // static census SP that needs to be preserved between timesteps
    typedef SP<Particle_Buffer<Particle<OS_Mesh> >::Census> SP_Census;
    static SP_Census host_census;

public:
  // constructor
    AMR_Interface(const Arguments &arg);
    
  // accessor functions

  // accessors for AMR_Builder
    const double* get_node_coord() const { return arguments.node_coord; }
    const int* get_layout() const { return arguments.layout; }
    int get_num_cells() const { return arguments.num_cells; }
    int get_num_global_cells() const { return arguments.global_num_cells; }

  // accessors for Opacity<MT>
    vector<double> get_density() const;
    vector<double> get_kappa() const;
    vector<double> get_kappa_thomson() const;
    vector<double> get_specific_heat() const;
    vector<double> get_temperature() const;
    double get_implicit() const { return arguments.implicitness; }
    string get_analytic_opacity() const { return "opacity"; }
    string get_analytic_sp_heat() const { return "dedt"; }

  // accessors for Source_Init<MT>
    const vector<double>& get_evol_ext() const { return evol_ext; }
    const vector<double>& get_rad_source() const { return evol_ext; }
    double get_rad_s_tend() const { return arguments.rad_s_tend; }
    const vector<string>& get_ss_pos() const { return ss_pos; } 
    const vector<double>& get_ss_temp() const { return ss_temp; } 
    vector<double> get_rad_temp() const;
    vector<int> get_global_cells() const;
    double get_delta_t() const { return arguments.delta_t; }
    double get_elapsed_t() const { return arguments.elapsed_t; }
    int get_npmax() const { return arguments.npmax; }
    int get_npnom() const { return arguments.npnom; }
    double get_dnpdt() const { return arguments.dnpdt; }
    string get_ss_dist() const { return "none"; }
    int get_printf() const { return arguments.print_f; }
    int get_buffer() const { return arguments.buffer; }
    int get_seed() const { return arguments.seed; }
    int get_capacity() const { return 100000; }
    int get_cycle() const { return arguments.cycle; }

  // T^4 slope comes from the interface here, we may need to work a more
  // "uniform" method out for general cases and stand-alone
    vector<vector<double> > get_t4_slope() const;

  // static functions to access in-between timestep variables
    static SP_Census get_census() { return host_census; }
    static void set_census(const SP_Census &cen) { host_census = cen; }
};

CSPACE

#endif                          // __imc_AMR_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of imc/AMR_Interface.hh
//---------------------------------------------------------------------------//
