//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/MC_Test.hh
 * \author Thomas M. Evans
 * \date   Wed Apr 21 19:33:36 1999
 * \brief  Components that we will use throughout the MC test suite.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_test_MC_Test_hh__
#define __mc_test_MC_Test_hh__

#include "../OS_Mesh.hh"
#include "../Topology.hh"
#include "../RZWedge_Mesh.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <vector>

namespace rtt_mc_test
{

//===========================================================================//
// Parser (OS_Mesh)
//===========================================================================//

class Parser
{
  private:
    std::string mesh_file;

  public:
    // constructor
    Parser() : mesh_file("OS_Input") {/*...*/}
    Parser(const std::string &f) : mesh_file(f) {/*...*/}
    
    // public copy functions for Mesh
    std::string get_mesh_file() const { return mesh_file; }
};

//===========================================================================//
// Topology class tests
//===========================================================================//
// Full replication topology test --> uses base class reference to derived
// class type

bool topology_replication_test(rtt_dsxx::SP<rtt_mc::OS_Mesh> mesh, 
			       const rtt_mc::Topology &top);

//---------------------------------------------------------------------------//
// Full DD topology test --> uses base class reference to derived
// class type

bool topology_DD_test(rtt_dsxx::SP<rtt_mc::OS_Mesh> mesh,
		      const rtt_mc::Topology &top);

//===========================================================================//
// MAKE AN AMR RZWEDGE_MESH
//===========================================================================//
// 12 cell amr rzwedge mesh:
/*
                     0
            ___________________
            |     |     |  |  |
	    |     |     |11|12|
            |  7  |  8  |__|__| 
	    |     |     |  |  |  
	r   |     |     |9 |10| (where r = x/(sqrt(phi/sinphi)cos(phi/2))
	    |_____|_____|__|__|
   grad = 0 |     |  |  |     | 0
	    |     |4 |5 |     |
            |  1  |__|__|  6  | 
	    |     |  |  |     |  
	    |     |2 |3 |     |
	    |_____|__|__|_____|
                
                     0
                     z
*/
rtt_dsxx::SP<rtt_mc::RZWedge_Mesh> make_RZWedge_Mesh_AMR(double phi);

typedef std::vector<double>    sf_double;
typedef std::vector<sf_double> vf_double;

//===========================================================================//
// PARTICLE CLASS FOR TESTING PARTICLE BUFFER
//===========================================================================//

int get_particle_size(rtt_rng::Rnd_Control);

class Dummy_Particle
{
  private:
    // Particle data.
    double                       w;
    int                          cell;

    rtt_dsxx::SP<rtt_rng::Sprng> random;

  public:
    // Constructor.
    Dummy_Particle(int, double, rtt_rng::Sprng);

    // Unpacking constructor.
    Dummy_Particle(const std::vector<char> &);

    // Packing function.
    std::vector<char> pack() const;

    // >>> ACCESSORS

    double                get_wt()   const { return w; }
    int                   get_cell() const { return cell; }
    const rtt_rng::Sprng& get_rn()   const { return *random; }

    // >>> TRANSPORT
    void transport(int = 0, double = 0.0, int = 0);

    // >>> MODIFIERS
    void set_cell(int c) { cell = c; }
};

} // end namespace rtt_mc_test

#endif                          // __mc_test_MC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of mc/test/MC_Test.hh
//---------------------------------------------------------------------------//
