#include "LAMG_F90.hh"
#include "ds++/Assert.hh"
#include <iostream>

namespace rtt_LAMG
{

extern "C" 
void LAMG_SOLVER_CONSTRUCT(RTT_F90_INTEGER &ptr,
			   const RTT_F90_INTEGER *iopts,
			   const RTT_F90_INTEGER &niopts,
			   const RTT_F90_DOUBLE *dopts,
			   const RTT_F90_INTEGER &ndopts,
			   RTT_F90_INTEGER &ierr)
{
    ptr = 42;
    std::cout << "In stubbed LAMG_SOLVER_CONSTRUCT"
	      << ", ptr = " << ptr << "." << std::endl;
    ierr = 0;
}
    
extern "C"
void LAMG_SOLVER_DESTRUCT(const RTT_F90_INTEGER &ptr,
			  RTT_F90_INTEGER &ierr)
{
    std::cout << "In stubbed LAMG_SOLVER_DESTRUCT("
	      << ptr << ")."<< std::endl;
    ierr = ptr == 42 ? 0 : -1;
}

extern "C"
void LAMG_SOLVER_SOLVE(const RTT_F90_INTEGER &ptr,
		       RTT_F90_DOUBLE *x,
		       const RTT_F90_INTEGER &nrows,
		       const RTT_F90_INTEGER &nentries,
		       const RTT_F90_DOUBLE *b,
		       const RTT_F90_INTEGER *rowPointer,
		       const RTT_F90_INTEGER *colIndex,
		       const RTT_F90_DOUBLE *val,
		       RTT_F90_INTEGER &ierr)
{
    std::cout << "In stubbed LAMG_SOLVER_SOLVE("
	      << ptr << ")." << std::endl;
    ierr = ptr == 42 ? 0 : -1;
}

} // end namespace rtt_LAMG

