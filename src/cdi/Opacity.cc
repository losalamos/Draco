//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/Opacity.cc
 * \author Kelly Thompson
 * \date   Thu Jun 23 13:57:07 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Opacity.hh"
#include "GandolfWrapper.hh"

#include <string>

namespace rtt_cdi
{
//------------------------------------------------------------------------------
//                        Opacity Base Class
//------------------------------------------------------------------------------


    // no contents.


//------------------------------------------------------------------------------
//                        Gandolf Opacities
//------------------------------------------------------------------------------

// Opacity constructor
//----------------------------------------------------------------------
GandolfOpacity::GandolfOpacity( string _data_filename )
    : dataFilename ( _data_filename )
	  
    {
	cout << "In GandolfOpacity::GandolfOpacity()" << endl;
	
	// Assert that the opacity data file exists and that it can be read.
	// Let Gandolf do this for now.
	//
	// ofstream infile( gandolfFilename );
	// Insist(!infile,"Could not open Gandalf data file for reading.");
	
	cout << "  calling gmatids via tg wrapper" << endl;

	const int kmat = 2;
	int matids[kmat], nmat, ier;

	gmatids( dataFilename, matids, kmat, nmat, ier );
    }
 

   // Return the inerpolated Rosseland Gray Opacity for the specified 
   // temperature and density.
   //----------------------------------------------------------------------
 double GandolfOpacity::getGray( const double temp, const double density )
     {
	 double grayOpacity = 0.0;
	 
	 return grayOpacity;
     }




} // end namespace rtt_cdi


//---------------------------------------------------------------------------//
//                              end of Opacity.cc
//---------------------------------------------------------------------------//
