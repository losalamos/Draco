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
#include "ds++/SP.hh"
#include <vector>

namespace rtt_mc_test
{

//===========================================================================//
// FAILURE LIMIT
//===========================================================================//

bool fail(int line);

bool fail(int line, char *file);

#define FAILURE fail(__LINE__, __FILE__);

//===========================================================================//
// Parser (OS_Mesh)
//===========================================================================//

class Parser
{
  public:
    // constructor
    Parser() {/*...*/}
    
    // public copy functions for Mesh
    std::string get_mesh_file() const { return "OS_Input"; }
};

//===========================================================================//
// Parser (RZWedge_Mesh)
//===========================================================================//

class Parser_RZ
{
  public:
    // constructor
    Parser_RZ() {/*...*/}
    
    // public copy functions for Mesh
    std::string get_mesh_file() const { return "RZWedge_Input"; }
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
	r   |     |     |9 |10| (where x = r/(sqrt(phi/sinphi)cos(phi/2))
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

} // end namespace rtt_mc_test

#endif                          // __mc_test_MC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of mc/test/MC_Test.hh
//---------------------------------------------------------------------------//
