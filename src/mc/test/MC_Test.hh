//----------------------------------*-C++-*----------------------------------//
// MC_Test.hh
// Thomas M. Evans
// Wed Apr 21 19:33:36 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Some free functions that we will use throughout the MC test suite
//---------------------------------------------------------------------------//

#ifndef __mc_test_MC_Test_hh__
#define __mc_test_MC_Test_hh__

#include "ds++/SP.hh"
#include <vector>
#include <iostream>
#include <string>

namespace rtt_mc_test
{

//===========================================================================//
// FAILURE LIMIT
//===========================================================================//

bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    return false;
}

//===========================================================================//
// 2D Objects Interface
//===========================================================================//
// make an interface for a 6 cell mesh

class MC_Interface
{
    typedef std::vector<std::vector<double> > vec_double;
    typedef std::vector<std::string> vec_string;

  private:
    std::string coord;
    vec_double fine_edge;
    vec_string bnd;

  public:
    // constructor
    inline MC_Interface();
    
    // public copy functions for Mesh
    std::string get_coordinates() const {return coord;}
    vec_double get_fine_edge() const {return fine_edge;} 
    vec_string get_boundaries() const {return bnd;}
};

// constructor
MC_Interface::MC_Interface()
    : coord("xy"), fine_edge(2), bnd(4)
{
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
}

} // end namespace rtt_mc_test

#endif                          // __mc_test_MC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of mc/test/MC_Test.hh
//---------------------------------------------------------------------------//
