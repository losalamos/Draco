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

//===========================================================================//
// F90 Functional Interface to Rage 
//===========================================================================//
// direct functional interface to F90 Rage code

extern "C"
{
    extern void rage_imc_(int *, double *, int *, int *, int *, int *);
}

//===========================================================================//
// class AMR_Interface
//===========================================================================//

IMCSPACE

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
	
	const int num_cells;
	const int num_b_cells;

      // constructor
	Arguments(const double *, const int *, const int *, const int *, int, 
		  int);
    };

private:
  // data from Rage
    Arguments arguments;

public:
  // constructor
    AMR_Interface(const Arguments &arg);
    
  // accessor functions
    const double* get_node_coord() const { return arguments.node_coord; }
    const int* get_layout() const { return arguments.layout; }
    int get_num_cells() const { return arguments.num_cells; }
};

CSPACE

#endif                          // __imc_AMR_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of imc/AMR_Interface.hh
//---------------------------------------------------------------------------//
