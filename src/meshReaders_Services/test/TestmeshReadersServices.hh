//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders_Services/test/TestmeshReadersServices.hh
 * \author B.T. Adams
 * \date   Tue Mar 14 09:48:00 2000
 * \brief  Header file for the meshReaders_Services class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __test_TestmeshReadersServices_hh__
#define __test_TestmeshReadersServices_hh__

#include "RTT_Format_Reader/RTT_Mesh_Reader.hh"
#include "../Connect.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <map>

namespace rtt_meshReaders_Services_test
{

typedef rtt_RTT_Format_Reader::RTT_Mesh_Reader RTT_Mesh_Reader;

enum Meshes {DEFINED, AMR};

extern std::vector<std::string>                        bndry_flags;
extern std::vector<int>                                flag_numbs;
extern rtt_dsxx::SP<RTT_Mesh_Reader>                   mesh;
extern rtt_dsxx::SP<rtt_meshReaders_Services::Connect> connect;

// ACCESSORS
bool check_connect(rtt_dsxx::SP<RTT_Mesh_Reader> mesh_, 
		   rtt_dsxx::SP<rtt_meshReaders_Services::Connect> connect_, 
		   const Meshes & meshtype_);


template<class VECTYPE>
struct compareXYZ
{
    bool operator()(const VECTYPE & low_val, const VECTYPE & high_val) const
    {
        // require agreement to six significant figures for equality. Note that
        // Shawn's tet mesh only agrees to four significant digits.
        const double EPSILON = 1.0e-06;
	std::vector<double> epsilon(low_val.size());
	bool sorted = true;

	Insist(low_val.size() == high_val.size(),"Improper sort arguments!");
	int dim = low_val.size();
	for (int d = dim - 1; d >= 0; d--)
	{
	    if (low_val[d] != 0 && high_val[d] != 0)
	        epsilon[d] = EPSILON * ((std::fabs(low_val[d]) + 
					 std::fabs(high_val[d]))/2.);
	    else
	        epsilon[d] = EPSILON;
	    // this strange looking logical operator will sort x,y,(z) 
	    // coordinates with x varying fastest, followed by y, and lastly z.
	    if (high_val[d] < low_val[d] && 
		std::fabs(low_val[d] - high_val[d]) > epsilon[d] && 
		(d == dim-1 || 
		 (d == dim-2 && 
		  std::fabs(high_val[d+1]-low_val[d+1]) < epsilon[d+1]) ||
		 (std::fabs(high_val[d+1]-low_val[d+1]) < epsilon[d+1]
		  && std::fabs(high_val[d+2]-low_val[d+2])<epsilon[d+2])))
	    {
	        sorted = false;
		d = -1;
	    }
	}
	return sorted;
    }
};

} // end namespace rtt_meshReaders_Services_test

#endif                // _test_TestmeshReadersServices_hh__

//---------------------------------------------------------------------------//
//               end of meshReaders/test/TestmeshReadersServices.hh
//---------------------------------------------------------------------------//
