//----------------------------------*-C++-*----------------------------------//
// pcg_DB.hh
// Dave Nystrom
// Fri May  9 09:45:08 1997
//---------------------------------------------------------------------------//
// @> pcg descriptor block.
//---------------------------------------------------------------------------//

#ifndef __linalg_pcg_DB_hh__
#define __linalg_pcg_DB_hh__

#include "nml/Group.hh"
#include "nml/Items.hh"

#include "ds++/String.hh"

//===========================================================================//
// class pcg_DB - pcg descriptor block

// Holds namelist parameters for setting up and controlling pcg.
//===========================================================================//

namespace rtt_PCGDiffusionSolver
{

class pcg_DB {
  public:
// Some way to name this thing.

    rtt_dsxx::String name;

    int itmeth;
    int nout;
    int levout;
    int itsmax;
    int ntest;
    int iqside;
    int iuinit;
    int ns1;
    int ns2;
    int ickstg;
    int iuexac;
    int idot;
    int istats;
    double ctimer;
    double rtimer;
    double flopsr;
    double zeta;
    double alpha;

// Methods.

  public:
    pcg_DB( rtt_dsxx::String _name ) : name(_name) {}

    pcg_DB(int itmeth_in, int levout_in, int itsmax_in, double zeta_in,
	   double alpha_in)
	: itmeth(itmeth_in), nout(6), levout(levout_in), itsmax(itsmax_in),
	  ntest(-1), iqside(0), iuinit(1), ns1(10), ns2(10), ickstg(-1),
	  iuexac(0), idot(1), istats(0), ctimer(0.0), rtimer(0.0),
	  flopsr(0.0), zeta(zeta_in), alpha(alpha_in), name("pcg_DB")
    {
	// empty
    }

    rtt_dsxx::String Name() const { return name; }

    void setup_namelist( NML_Group& g );
};

} // end namespace rtt_PCGDiffusionSolver

#endif                          // __linalg_pcg_DB_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/pcg_DB.hh
//---------------------------------------------------------------------------//
