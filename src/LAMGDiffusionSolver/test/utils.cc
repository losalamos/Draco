//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/test/utils.cc
 * \author Randy M. Roberts
 * \date   Wed Oct 27 14:30:54 1999
 * \brief  utilities for testing the LAMGDiffusionSolver
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "utils.hh"
#include "Mesh_XYZFactory.hh"
#include "../Release.hh"
#include "nml/Group.hh"

#include <iostream>

namespace rtt_LAMGDiffusionSolver_test
{

void version(const std::string &progname)
{
    std::cout << "LAMGDiffusionSolver: version "
	      << rtt_LAMGDiffusionSolver::release() << std::endl;
}

Mesh_XYZFactory getMTFactory(const std::string &filename)
{
    NML_Group g( "test" );
    Mesh_XYZFactory::Mesh_DB mdb;
    mdb.setup_namelist( g );
    g.readgroup( filename.c_str() );
    
    std::cout << "Factory created in Tester" << std::endl;

    return Mesh_XYZFactory(mdb);
}

 
} // end namespace rtt_LAMGDiffusionSolver_test


//---------------------------------------------------------------------------//
//                              end of utils.cc
//---------------------------------------------------------------------------//
