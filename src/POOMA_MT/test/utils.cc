//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   POOMA_MT/test/utils.cc
 * \author Randy M. Roberts
 * \date   Wed Oct 27 14:30:54 1999
 * \brief  utilities for testing the POOMA_MT
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "utils.hh"

#include "meshTest/Release.hh"
#include "../Release.hh"

#include <fstream>

namespace
{
inline int Dimension() { return 3; }
}

namespace rtt_POOMA_MT_test
{

void version(const std::string &progname)
{
    cout << "POOMA_MT: version " << rtt_POOMA_MT::release() << endl;
    cout << "meshTest: version " << rtt_meshTest::release() << endl;
}

PoomaMesh_XYZFactory getMTFactory(const std::string &filename)
{
    std::ifstream ifs(filename.c_str());
    
    std::vector<int> numCells(Dimension());
    for (int i=0; i<Dimension(); i++)
	ifs >> numCells[i];

    std::vector<std::vector<double> > cellWidth(Dimension());
    for (int i=0; i<Dimension(); i++)
    {
	cellWidth[i] = std::vector<double>(numCells[i]);
	for (int j=0; j<numCells[i]; j++)
	    ifs >> cellWidth[i][j];
    }

    return PoomaMesh_XYZFactory(numCells, cellWidth);
}

 
} // end namespace rtt_POOMA_MT_test


//---------------------------------------------------------------------------//
//                              end of utils.cc
//---------------------------------------------------------------------------//
