//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/GandolfWrapper.cc
 * \author Kelly Thompson
 * \date   Thu Jun 29 15:31:56 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfWrapper.hh"

#include <string>

namespace rtt_cdi {

    //----------------------------------------//
    //                gmatids                 //
    //----------------------------------------//
    
    void gmatids( const std::string &fname , int matids[], 
		  const int kmat, int &nmat, int &ier ) 
	{
	    int len = fname.length();
	    // assert( len < 80 )

	    // copy filename into a char array;
	    // also remove constness.
	    char cfname[80];
	    int i = 0;
	    for ( ; i < len; ++i )
		cfname[i] =fname[i];
	    for ( ; i < 80; ++i )
		cfname[i] = ' ';

	    // create "long int" versions of variables.
	    long int li_nmat = static_cast<long int>(nmat);
	    long int li_ier  = static_cast<long int>(ier);
	    // also remove constness from kmat.
	    long int li_kmat = static_cast<long int>(kmat); 

	    // we don't know the value of kmat until runtime so we
	    // must dynamically allocate li_matids.
	    long int *li_matids = new long int ( kmat );
	    for ( int i=0; i<kmat; ++i )
		li_matids[i] = static_cast<long int>(matids[i]);

	    // call the Gandolf library function
	    extc_gmatids( cfname, li_matids, li_kmat, li_nmat, li_ier );

	    // update the function arguments from their "long int"
	    // counterparts.
	    for ( int i=0; i<kmat; ++i )
		matids[i] = static_cast<int>(li_matids[i]);
	    ier  = static_cast<int>(li_ier);
	    nmat = static_cast<int>(li_nmat);
	    // we don't update kmat since it is const.
	    delete li_matids;

	    return;
	}


} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
//                         end of cdi/GandolfWrapper.cc
//---------------------------------------------------------------------------//
