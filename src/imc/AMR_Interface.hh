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
#include "ds++/SP.hh"

//===========================================================================//
// F90 Functional Interface to Rage 
//===========================================================================//
// direct functional interface to F90 Rage code

extern "C"
{
    extern void rage_imc_(int *, double *, int *, int *, int *, int *);
}

//---------------------------------------------------------------------------//
// census object declarations

IMCSPACE
GLOBALSPACE

extern 
dsxx::SP<IMC::Particle_Buffer<IMC::Particle<IMC::OS_Mesh> >::Census>
host_census; 

CSPACE

//===========================================================================//
// class AMR_Interface
//===========================================================================//

class AMR_Interface
{
public:
  // structure for passing arguments to the interface
    struct Arguments
    {
      // data determining the problem layout
	const double *node_coord;
	const int *layout;
	const int *b_proc;
	const int *b_cell;
	const double *dedt;
	const double *rho;
	const double *opacity_abs;
	const double *tev;
	const double *rev;
	
	const int num_cells;
	const int num_b_cells;	
	const double implicitness;
	const double delta_t;
	const double dnpdt;
	const int npnom;
	const int npmax;
	const double rad_s_tend;
	const int seed;
	const int buffer;
	const int print_f;

      // constructor
	Arguments(const double *, const int *, const int *, const int *,
		  const double *, const double *, const double *, const
		  double *, const double *, int, int, double, double, double, 
		  int, int, double, int, int, int); 
    };

private:
  // data from Rage
    Arguments arguments;

  // blank vector argument to be filled in at later date when the proper
  // features are added
    vector<double> evol_ext;
    vector<string> ss_pos;
    vector<double> ss_temp;

public:
  // constructor
    AMR_Interface(const Arguments &arg);
    
  // accessor functions

  // accessors for AMR_Builder
    const double* get_node_coord() const { return arguments.node_coord; }
    const int* get_layout() const { return arguments.layout; }
    int get_num_cells() const { return arguments.num_cells; }

  // accessors for Opacity<MT>
    vector<double> get_density() const;
    vector<double> get_kappa() const;
    vector<double> get_specific_heat() const;
    vector<double> get_temperature() const;
    double get_implicit() const { return arguments.implicitness; }

  // accessors for Source_Init<MT>
    const vector<double>& get_evol_ext() const { return evol_ext; }
    const vector<double>& get_rad_source() const { return evol_ext; }
    double get_rad_s_tend() const { return arguments.rad_s_tend; }
    const vector<string>& get_ss_pos() const { return ss_pos; } 
    const vector<double>& get_ss_temp() const { return ss_temp; } 
    vector<double> get_rad_temp() const;
    double get_delta_t() const { return arguments.delta_t; }
    int get_npmax() const { return arguments.npmax; }
    int get_npnom() const { return arguments.npnom; }
    double get_dnpdt() const { return arguments.dnpdt; }
    string get_ss_dist() const { return "none"; }
    string get_analytic_opacity() const { return "opacity"; }
    string get_analytic_sp_heat() const { return "dedt"; }
    int get_printf() const { return arguments.print_f; }
    int get_buffer() const { return arguments.buffer; }
    int get_seed() const { return arguments.seed; }
    int get_capacity() const { return 0; }
};

CSPACE

#endif                          // __imc_AMR_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of imc/AMR_Interface.hh
//---------------------------------------------------------------------------//
