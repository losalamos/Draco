//----------------------------------*-C++-*----------------------------------//
// utils.t.cc
// Randy M. Roberts
// Wed Nov  4 13:47:25 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/utils.hh"
#include "c4/global.hh"

namespace rtt_3T_testP13T
{

 template<class MT>
 bool isContinuous(const typename MT::fcdsf &rhs,
		   const SP<MT> &spMesh)
 {
     Require(*spMesh == rhs.get_Mesh());
     
     typename MT::fcdsf swapped(spMesh);

     // We test continuity by swapping faces across cells and compare
     // as negative values.
     
     MT::swap(swapped, rhs);

     // We don't want to have to worry about boundaries, so I strip them
     // off of the original field, invert them, and assign them to the
     // swapped field to simulate continuitiy at the boundary.
     
     typename MT::bssf boundary(spMesh);
     MT::gather(boundary, rhs, MT::OpAssign());
     boundary *= -1.0;
     MT::gather(swapped, boundary, MT::OpAssign());

     int nodeContinuous = 1;
     
     for (MT::fcdsf::const_iterator rit = rhs.begin(),
	      sit = swapped.begin();
	  rit != rhs.end(); rit++, sit++)
     {
	 if (*rit != -*sit)
	 {
	     nodeContinuous = 0;
	     break;
	 }
     }

     int numNodesContinuous = nodeContinuous;

     C4::gsum(numNodesContinuous);
     return numNodesContinuous == C4::nodes();
 }
 
 template<class FT>
 double sum(const FT &rhs)
 {
     double results = 0.0;
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	 results += *it;

     C4::gsum(results);
     return results;
 }

 template<class FT>
 typename FT::value_type min(const FT &rhs)
 {
     FT::value_type results = *rhs.begin();
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	 if (*it < results)
	     results = *it;

     C4::gmin(results);
     return results;
 }

 template<class FT>
 typename FT::value_type max(const FT &rhs)
 {
     FT::value_type results = *rhs.begin();
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	 if (*it > results)
	     results = *it;

     C4::gmax(results);
     return results;
 }

 template<class FT, class OP>
 typename FT::value_type min(const FT &rhs, OP &op)
 {
     FT::value_type results = *rhs.begin();
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	 if (op(*it) < results)
	     results = *it;

     C4::gmin(results);
     
     return results;
 }

 template<class FT, class OP>
 typename FT::value_type max(const FT &rhs, OP &op)
 {
     FT::value_type results = *rhs.begin();
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
	 if (op(*it) > results)
	     results = *it;

     C4::gmax(results);
     
     return results;
 }

} // end namespace rtt_3T_testP13T

//---------------------------------------------------------------------------//
//                              end of utils.t.cc
//---------------------------------------------------------------------------//
