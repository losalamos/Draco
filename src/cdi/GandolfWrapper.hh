//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/GandolfWrapper.hh
 * \author Kelly Thompson
 * \date   Thu Jun 29 15:31:56 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_GandolfWrapper_hh__
#define __cdi_GandolfWrapper_hh__

#include <string>

// function prototypes
//---------------------------------------------------------------------------//

namespace rtt_cdi {
 
    void gmatids( const std::string &fname , int matids[], 
		  const int kmat,int &nmat, int &ier );

// I need to create a "long int" version of the calling routine.

//     void gmatids( const std::string &fname , long int matids[], 
// 		  const long int kmat, long int &nmat, long int &ier );

} // end namespace rtt_cdi



// Handle machine specific FORTRAN name linkage.
//---------------------------------------------------------------------------//

#if defined(sun) || defined(__sun) || defined(__sgi) || defined(__linux)
    
#define extc_gmatids gmatids_
    
#endif
    


// Function prototypes for Gandolf F77 subroutines.
//---------------------------------------------------------------------------//

// The Gandolf library was compiled with -i8 so we must use "long int" 
// values.

extern "C" {
    void extc_gmatids( char *cfname, long int *matids, long int &ckmat,
		       long int &nmat, long int &ier );

    // This call may be used for calling a test routine ( with normal
    // ints instead of longs ).
    // ----------------------------------------------------------------
//     void extc_gmatids( char *cfname, int *matids, int &ckmat,
// 		       int &nmat, int &ier );
}



#endif                          // __cdi_GandolfWrapper_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/GandolfWrapper.hh
//---------------------------------------------------------------------------//
