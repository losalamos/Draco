//----------------------------------*-C++-*----------------------------------//
// GmvDump.hh
// Randy M. Roberts
// Wed Oct 14 10:19:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_testP13T_GmvDump_hh__
#define __3T_testP13T_GmvDump_hh__

#include "ds++/SP.hh"
#include <string>

namespace rtt_3T_testP13T
{

 using dsxx::SP;
 
 //===========================================================================//
 // class GmvDump - 
 //
 // Date created :
 // Purpose      :
 //
 // revision history:
 // -----------------
 // 0) original
 // 
 //===========================================================================//

 template<class MT>
 class GmvDump
 {

     // NESTED CLASSES AND TYPEDEFS

   private:
     
     struct VertId
     {
	 int nx;
	 int ny;
	 int nz;
	 VertId(int nx_, int ny_, int nz_) : nx(nx_+1), ny(ny_+1), nz(nz_+1) { }
	 int operator()(int i, int j, int k)
	 {
	     return k*nx*ny + j*nx + i + 1;
	 }
     };
     
     // DATA
    
     std::string fname;
     const SP<MT> &spMesh;

     int nx;
     int ny;
     int nz;
     int nzp;
     int zoff;

     bool variablePrinted;

     int cycle;
     double time;
     
   public:

     // CREATORS
    
     GmvDump(const std::string &fname_, const SP<MT> &spMesh_,
	     int cycle_, double time_);

     ~GmvDump();
     
     // MANIPULATORS

     void dump(const typename MT::ccsf &var, const std::string &name);
     
     // ACCESSORS

   private:
    
     // IMPLEMENTATION
 };

} // end namespace rtt_3T_testP13T

#endif                          // __3T_testP13T_GmvDump_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/GmvDump.hh
//---------------------------------------------------------------------------//
